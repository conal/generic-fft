{-# LANGUAGE ScopedTypeVariables, ConstraintKinds #-}
{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
{-# LANGUAGE UndecidableInstances #-}      -- See below
{-# OPTIONS_GHC -Wall #-}

----------------------------------------------------------------------
-- |
-- Module      :  Generic.FFT
-- Copyright   :  (c) 2013 Tabula, Inc.
-- 
-- Maintainer  :  conal@tabula.com
-- Stability   :  experimental
-- 
-- FFT via functor composition. Apply to depth-typed perfect binary leaf trees.
----------------------------------------------------------------------

module Generic.FFT (dft,HasFFT(..)) where

import Prelude hiding (sum)

import Data.Monoid (Monoid(..),Sum(..),(<>))
import Data.Functor ((<$>))
import Data.Foldable (Foldable,sum)
import Data.Traversable (Traversable(..), mapAccumL)
import Control.Applicative (Applicative(..),liftA2)
import Control.Arrow ((***))

import Data.Tuple (swap)
import Data.Complex (Complex(..))

import Data.UniformPair

import TypeUnary.Nat

import Data.FTree.BottomUp

-- TODO: Explore top-down trees as well.

{--------------------------------------------------------------------
    Misc
--------------------------------------------------------------------}

type Unop a = a -> a

transpose :: (Traversable g, Applicative f) => g (f a) -> f (g a)
transpose = sequenceA

type C = Complex Double

scanL :: (Traversable f, Monoid a) => (a, f a) -> (f a, a)
scanL (a0,as) = swap (mapAccumL h a0 as)
 where
   h a a' = (b,b) where b = a <> a'

-- TODO: Replace scanL with my efficient parallel version.

-- Prefix (left) sums
sumsL :: (Traversable f, Num a) => (a, f a) -> (f a, a)
sumsL = (fmap getSum *** getSum) . scanL . (Sum *** fmap Sum)

-- Yield a structure counting from 0 to size-1, together with size
counts :: forall f a. (TA f, Num a) => (f a, a)
counts = sumsL (0, pure 1)

-- Cross product of structures
cross :: (Functor g, Functor f) => g a -> f b -> g (f (a,b))
as `cross` bs = fmap (\ a -> fmap (\ b -> (a,b)) bs) as

-- All products of numbers from each structure.
products :: (Functor g, Functor f, Num n) => g n -> f n -> g (f n)
products = (fmap.fmap.fmap.fmap) (uncurry (*)) cross

-- as `products` bs = (fmap.fmap) (uncurry (*)) (as `cross` bs)
-- as `products` bs = fmap (\ a -> fmap (\ b -> a*b) bs) as

-- Dot product of structures
dot :: (Applicative f, Foldable f, Num a) => f a -> f a -> a
u `dot` v = sum (liftA2 (*) u v)

i2pi :: C
i2pi = 0 :+ 2*pi

-- Principle nth root of unity with negated angle
uroot :: Int -> C
uroot n = exp (- i2pi / fromIntegral n)

{--------------------------------------------------------------------
    Unoptimized discrete Fourier transform (DFT)
--------------------------------------------------------------------}

-- Discrete Fourier transform:
-- $X_k = \sum_{n=0}^{N-1} x_n e^{-i 2\pi k \frac{n}{N}}$ for $k = 0,\ldots,N$.
dft :: TA f => Unop (f C)
dft xs = (xs `dot`) <$> rootses

-- Powers of 'uroot' needed in the DFT:
-- $e^{\frac{-i 2\pi k n}{N}}$ for $k,n = 0,\ldots,N$:
rootses :: forall f. TA f => f (f C)
rootses = rootCross tot indices indices
 where
   indices :: f Int
   (indices,tot) = counts

rootCross :: (Functor g, Functor f, Integral n) =>
             Int -> g n -> f n -> g (f C)
rootCross tot is js = (fmap.fmap) (uroot tot ^) (is `products` js)

{--------------------------------------------------------------------
    FFT
--------------------------------------------------------------------}

-- | FFT computation, parametrized by structure
class HasFFT f where
  fft :: Unop (f C)

-- Constraint shorthands
type TA  f = (Traversable f, Applicative f)
type TAH f = (TA f, HasFFT f)

instance HasFFT Pair where
  fft (a :# b) = a+b :# a-b

instance (TAH f, IsNat n) => HasFFT (T f n) where
  fft = inT id fftC

--     Variable s `f, f' occur more often than in the instance head
--       in the constraint: TAH f
--     (Use -XUndecidableInstances to permit this)
--     In the instance declaration for `HasFFT (T f n)'
-- 
-- This warning vanishes when we spell out TAH. Hm.

-- FFTs after transposition
fftsT :: (Applicative f, Traversable g, HasFFT g) => g (f C) -> f (g C)
fftsT = fmap fft . transpose

--   transpose :: g (f C) -> f (g C)
--   fmap fft  :: f (g C) -> f (g C)

-- FFT of composed functors
fftC :: (TAH f, TAH g) => Unop (g (f C))
fftC = fftsT . twiddle . fftsT

--   fftsT   :: g (f C) -> f (g C)
--   twiddle :: f (g C) -> f (g C)
--   fftsT   :: f (g C) -> g (f C)

-- Multiply by twiddle factors
twiddle :: forall f g. (TA f, TA g) => Unop (g (f C))
twiddle = (liftA2.liftA2) (*) rootses'
 where
   rootses' = rootCross (gTot*fTot) gIndices fIndices
    where
      fIndices :: f Int
      (fIndices,fTot) = counts
      gIndices :: g Int
      (gIndices,gTot) = counts
