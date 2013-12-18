{-# LANGUAGE ScopedTypeVariables, ConstraintKinds, GADTs #-}
{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
{-# LANGUAGE UndecidableInstances #-}      -- See below

{-# OPTIONS_GHC -Wall #-}
{-# OPTIONS_GHC -fno-warn-unused-imports #-}
{-# OPTIONS_GHC -fno-warn-unused-binds   #-}

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

import qualified Data.FTree.BottomUp as B
import qualified Data.FTree.TopDown  as T

{--------------------------------------------------------------------
    Misc
--------------------------------------------------------------------}

type Unop a = a -> a

transpose :: (Traversable g, Applicative f) => g (f a) -> f (g a)
transpose = sequenceA

inTranspose :: (Traversable f, Traversable k, Applicative g, Applicative h) =>
               (g (f a) -> k (h b)) -> (f (g a) -> h (k b))
inTranspose = transpose --> transpose

infixr 1 -->
-- | Add pre- and post processing
(-->) :: (a' -> a) -> (b -> b') -> ((a -> b) -> (a' -> b'))
(f --> h) g = h . g . f

type C = Complex Double

scanL :: (Traversable f, Monoid a) => (a, f a) -> (f a, a)
scanL (a0,as) = swap (mapAccumL h a0 as)
 where
   h a a' = (a <> a',a)

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

-- Dot product of structures. Assumes a trie-like f, so that the Applicative
-- instance combines corresponding elements. Perhaps replace Applicative with a
-- more suitable constraint.
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
rootCross tot = (fmap.fmap.fmap.fmap) (uroot tot ^) products

-- rootCross tot is js = (fmap.fmap) (uroot tot ^) (is `products` js)

{--------------------------------------------------------------------
    FFT
--------------------------------------------------------------------}

-- | FFT computation, parametrized by structure
class HasFFT f f' | f -> f' where
  fft :: f C -> f' C

-- Constraint shorthands
type TA  f = (Traversable f, Applicative f)
type TAH f f' = (TA f, TA f', HasFFT f f')

-- Binary butterfly.
instance HasFFT Pair Pair where
  fft (a :# b) = a+b :# a-b

-- Decimation in time (DIT)
instance (TAH f f', IsNat n) => HasFFT (B.T f n) (T.T f' n) where
  fft (B.L a) = T.L a
  fft (B.B t) = T.B (fftC t)

-- Decimation in frequency (DIF)? I'm unsure.
instance (TAH f f', IsNat n) => HasFFT (T.T f n) (B.T f' n) where
  fft (T.L a) = B.L a
  fft (T.B t) = B.B (fftC t)

--     Variable s `f, f' occur more often than in the instance head
--       in the constraint: TAH f
--     (Use -XUndecidableInstances to permit this)
--     In the instance declaration for `HasFFT (T f n)'
-- 
-- This warning vanishes when we spell out TAH. Hm.

-- I'd prefer terser definitions like `fft = T.inT' T.l (T.B . fftC)`, but I
-- haven't found a type for `inT'` that GHC likes.

fftsT :: (Applicative f, Traversable g, HasFFT g g') => g (f C) -> f (g' C)
fftsT = fmap fft . transpose

-- FFT of composed functors
fftC :: (TAH f f', TAH g g') => g (f C) -> f' (g' C)
fftC = transpose . fftsT . twiddle . fftsT 

{-

--   fftsT     :: g  (f  C) -> f  (g' C)
--   twiddle   :: f  (g' C) -> f  (g' C)
--   fftsT     :: f  (g' C) -> g' (f' C)
--   transpose :: g' (f' C) -> f' (g' C)

-}

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

{--------------------------------------------------------------------
    Experimental variation
--------------------------------------------------------------------}

-- This version is a better fit with
-- <https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#General_factorizations>

-- FFT of composed functors
fftC' :: (TAH f f', TAH g g') => g (f C) -> f' (g' C)
fftC' = fftsT' . transpose . twiddle . fftsT'
 where
   fftsT' :: (TA h, TAH k k') => k (h C) -> k' (h C)
   fftsT' = (inTranspose.fmap) fft

-- Hm. This definition differs from the previous one, since `twiddle` and
-- `transpose` got swapped. I doubt they're equivalent.

-- Types in the fftC' definition (right to left):
--
--   fftsT'    :: g  (f  C) -> g' (f  C)
--   twiddle   :: g' (f  C) -> g' (f  C)
--   transpose :: g' (f  C) -> f  (g' C)
--   fftsT'    :: f  (g' C) -> g' (g' C)


{--------------------------------------------------------------------
    Tests
--------------------------------------------------------------------}

-- {1, 0, 0, 0, ...} <=DFT=> {1, 1, 1, 1, ...}\ 
-- {1, -1, 1, -1, 1, -1, ...} <=DFT=> {0, 0, ... , N, 0, 0, ...} (where the 'N' occurs at the N/2 position in the output vector.

_p1 :: Pair C
_p1 = 1 :# 0

_t1 :: B.T Pair N3 C
_t1 = B.B (B.B (B.B (B.L (((1 :# 0) :# (0 :# 0)) :# ((0 :# 0) :# (0 :# 0))))))

_t2 :: B.T Pair N4 C
_t2 = pure 1

_t3 :: B.T Pair N3 C
_t3 = B.B (pure (1 :# -1))
