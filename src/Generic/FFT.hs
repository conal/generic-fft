{-# LANGUAGE TypeOperators, GADTs #-}
{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
{-# LANGUAGE TypeSynonymInstances, FlexibleInstances, ScopedTypeVariables #-}
{-# LANGUAGE NoMonomorphismRestriction #-} -- TEMP
{-# LANGUAGE ConstraintKinds #-}           -- experimental
{-# LANGUAGE UndecidableInstances #-}      -- see below
{-# OPTIONS_GHC -Wall #-}

-- {-# OPTIONS_GHC -fno-warn-unused-imports #-} -- TEMP
-- {-# OPTIONS_GHC -fno-warn-unused-binds   #-} -- TEMP

----------------------------------------------------------------------
-- |
-- Module      :  Generic.FFT
-- Copyright   :  (c) 2013 Tabula, Inc.
-- 
-- Maintainer  :  conal@tabula.com
-- Stability   :  experimental
-- 
-- FFT on depth-typed perfect binary leaf trees
----------------------------------------------------------------------

module Generic.FFT where

-- TODO: explicit exports

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

-- Constraint shorthands
type TA  f = (Traversable f, Applicative f)
type TAH f = (TA f, HasFFT f)

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

counts :: forall f a. (TA f, Num a) => (f a, a)
counts = sumsL (0, pure 1)

products :: (Functor g, Functor f, Num n) => g n -> f n -> g (f n)
as `products` bs = (\ a -> (a *) <$> bs) <$> as

dot :: (Applicative f, Foldable f, Num a) => f a -> f a -> a
u `dot` v = sum (liftA2 (*) u v)

i2pi :: C
i2pi = 0 :+ 2*pi

uroot :: Int -> C
uroot n = exp (- i2pi / fromIntegral n)

{--------------------------------------------------------------------
    Unoptimized discrete Fourier transform (DFT)
--------------------------------------------------------------------}

-- $X_k = \sum_{n=0}^{N-1} x_n e^{-i 2\pi k \frac{n}{N}}$ for $k = 0,\ldots,N$.

dft :: TA f => Unop (f C)
dft xs = (xs `dot`) <$> rootses

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

instance HasFFT Pair where
  fft (a :# b) = a+b :# a-b

instance (TAH f, IsNat n) => HasFFT (T f n) where
  fft = inT id fftC

--     Variable s `f, f' occur more often than in the instance head
--       in the constraint: TAH f
--     (Use -XUndecidableInstances to permit this)
--     In the instance declaration for `HasFFT (T f n)'

ffts' :: (Applicative f, Traversable g, HasFFT g) => g (f C) -> f (g C)
ffts' = fmap fft . transpose

--   transpose :: g (f C) -> f (g C)
--   fmap fft  :: f (g C) -> f (g C)

fftC :: (TAH f, TAH g) => Unop (g (f C))
fftC = ffts' . twiddle . ffts'

--   ffts'   :: g (f C) -> f (g C)
--   twiddle :: f (g C) -> f (g C)
--   ffts'   :: f (g C) -> g (f C)

twiddle :: forall f g. (TA f, TA g) => Unop (g (f C))
twiddle = (liftA2.liftA2) (*) rootses'
 where
   rootses' = rootCross (gTot*fTot) gIndices fIndices
    where
      fIndices :: f Int
      (fIndices,fTot) = counts
      gIndices :: g Int
      (gIndices,gTot) = counts
