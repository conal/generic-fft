{-# LANGUAGE TypeOperators, GADTs #-}
{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}
{-# LANGUAGE TypeSynonymInstances, FlexibleInstances, ScopedTypeVariables #-}
{-# LANGUAGE NoMonomorphismRestriction #-} -- TEMP
{-# LANGUAGE ConstraintKinds #-}  -- experimental
{-# OPTIONS_GHC -Wall #-}

{-# OPTIONS_GHC -fno-warn-unused-imports #-} -- TEMP
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

import Data.Monoid (Monoid(..),Product(..),Sum(..))
import Data.Functor ((<$>))
import Data.Foldable (Foldable,sum)
import Data.Traversable (Traversable(..), mapAccumL)
import Control.Applicative (Applicative(..),liftA2)
import Control.Arrow ((***))

import Data.Complex (Complex(..))

import Data.UniformPair

import TypeUnary.Nat

import Data.FTree.BottomUp
-- import Data.FTree.TopDown

{--------------------------------------------------------------------
    Unoptimized discrete Fourier transform (DFT)
--------------------------------------------------------------------}

-- $X_k = \sum_{n=0}^{N-1} x_n e^{-i 2\pi k \frac{n}{N}}$ for $k = 0,\ldots,N$.

dft :: forall f . (Traversable f, Applicative f) => Unop (f C)
dft xs = (xs `dot`) <$> rootses

-- $e^{\frac{-i 2\pi k n}{N}}$:

rootses :: forall f. (Traversable f, Applicative f) => f (f C)
-- rootses = (fmap.fmap) ((uroot tot ^) . uncurry (*)) (indices `cross` indices)
rootses = rootCross tot indices indices
 where
   indices :: f Int
   (indices,tot) = counts

-- Experimental generalization for twiddle

rootCross :: (Functor g, Functor f, Integral n) =>
             Int -> g n -> f n -> g (f C)
rootCross tot is js = (fmap.fmap) ((uroot tot ^) . uncurry (*)) (is `cross` js)

{--------------------------------------------------------------------
    FFT
--------------------------------------------------------------------}

-- | FFT computation, parametrized by structure
class HasFFT f where
  fft :: Unop (f C)

instance HasFFT Pair where
  fft (a :# b) = a+b :# a-b

instance (Traversable f, Applicative f, HasFFT f, IsNat n)
      => HasFFT (T f n) where
  fft = inT id fftC

ffts' :: (Applicative f, Traversable g, HasFFT g) => g (f C) -> f (g C)
ffts' = fmap fft . transpose

--   transpose :: g (f C) -> f (g  C)
--   fmap fft  :: f (g C) -> f (g' C)

fftC :: ( Traversable f, Applicative f, HasFFT f
        , Traversable g, Applicative g, HasFFT g ) =>
        Unop (g (f C))
fftC = ffts' . twiddle . ffts'

-- Types:
-- 
--   ffts'   :: g (f C) -> f (g C)
--   twiddle :: f (g C) -> f (g C)
--   ffts'   :: f (g C) -> g (f C)

twiddle :: (Traversable f, Applicative f, Traversable g, Applicative g) =>
           Unop (g (f C))
twiddle = (liftA2.liftA2) (*) rootses'

rootses' :: forall g f. (Traversable f, Applicative f, Traversable g, Applicative g) => g (f C)
rootses' = rootCross (gTot*fTot) gIndices fIndices
 where
   fIndices :: f Int
   (fIndices,fTot) = counts
   gIndices :: g Int
   (gIndices,gTot) = counts

-- TODO: Factor out ((uroot tot ^) . uncurry (*)) from rootses and twiddle.
-- Does twiddle generalize rootses?


spread :: forall f. (Applicative f, Traversable f) => f C
spread = getProduct <$> fst (scanL (Product 1, pure (Product delta)))
 where
   delta = exp (i2pi / fromIntegral n)
   n     = sum (pure 1 :: f Int)

-- TODO: Package up parallel scan and use here.


{--------------------------------------------------------------------
    Misc
--------------------------------------------------------------------}

type Unop a = a -> a

type TA  f = (Traversable f, Applicative f)
type TAH f = (TA f, HasFFT f)

transpose :: (Traversable g, Applicative f) => g (f a) -> f (g a)
transpose = sequenceA

{-
inTranspose :: (Traversable f, Traversable k, Applicative g, Applicative h) =>
               (g (f a) -> k (h b)) -> (f (g a) -> h (k b))
inTranspose = sequenceA --> sequenceA

infixr 1 -->
-- | Add pre- and post processing
(-->) :: (a' -> a) -> (b -> b') -> ((a -> b) -> (a' -> b'))
(f --> h) g = h . g . f
-}

type C = Complex Double

scanL :: (Traversable f, Monoid a) => (a, f a) -> (f a, a)
scanL = undefined -- for now

-- Prefix (left) sums
sumsL :: (Traversable f, Num a) => (a, f a) -> (f a, a)
sumsL = (fmap getSum *** getSum) . scanL . (Sum *** fmap Sum)

counts :: forall f a. (Traversable f, Applicative f, Num a) => (f a, a)
counts = sumsL (0, pure 1)

cross :: (Functor g, Functor f) => g a -> f b -> g (f (a , b))
as `cross` bs = fmap (\ a -> fmap (\ b -> (a,b)) bs) as

dot :: (Applicative f, Foldable f, Num a) => f a -> f a -> a
u `dot` v = sum (liftA2 (*) u v)

i2pi :: C
i2pi = 0 :+ 2*pi

uroot :: Int -> C
uroot n = exp (- i2pi / fromIntegral n)

