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

-- dift :: forall f t. (Traversable f, RealFloat t) => Unop (f (Complex t))
-- dift = ...
--  where
--    twiddles = 
--    (indices,tot) = counts

{--------------------------------------------------------------------
    FFT
--------------------------------------------------------------------}

-- | FFT computation, parametrized by structure
class HasFFT f where
  fft :: RealFloat t => Unop (f (Complex t))

instance HasFFT Pair where
  fft (a :# b) = a+b :# a-b

instance (Traversable f, Applicative f, HasFFT f, IsNat n)
      => HasFFT (T f n) where
  fft = inT id fftC

fftC :: ( Traversable f, Traversable g, Applicative f, Applicative g
        , HasFFT f, HasFFT g, RealFloat t ) =>
        Unop (g (f (Complex t)))
fftC = fmap fft . transpose . twiddle . fmap fft . transpose

-- Types:
-- 
--   transpose :: g  (f  c) -> f  (g  c)
--   fmap fft  :: f  (g  c) -> f  (g' c)
--   twiddle   :: f  (g' c) -> f  (g' c)
--   transpose :: f  (g' c) -> g' (f  c)
--   fmap fft  :: g' (f  c) -> g' (f' c)

twiddle :: Unop (g (f (Complex t)))
twiddle = id

spread :: forall f t. (Applicative f, Traversable f, RealFloat t) => f (Complex t)
spread = getProduct <$> fst (scanL (Product 1, pure (Product delta)))
 where
   delta :: Complex t
   delta = exp (i2pi / fromIntegral n)
   n :: Int
   n = sum (pure 1 :: f Int)

i2pi :: RealFloat t => Complex t
i2pi = 0 :+ 2*pi

-- TODO: Package up parallel scan and use here.


{--------------------------------------------------------------------
    Misc
--------------------------------------------------------------------}

type Unop a = a -> a

transpose :: (Traversable g, Applicative f) => g (f a) -> f (g a)
transpose = sequenceA

inTranspose :: (Traversable f, Traversable k, Applicative g, Applicative h) =>
               (g (f a) -> k (h b)) -> (f (g a) -> h (k b))
inTranspose = sequenceA --> sequenceA


infixr 1 -->
-- | Add pre- and post processing
(-->) :: (a' -> a) -> (b -> b') -> ((a -> b) -> (a' -> b'))
(f --> h) g = h . g . f


{--------------------------------------------------------------------
    Experiments
--------------------------------------------------------------------}

scanL :: (Traversable f, Monoid a) => (a, f a) -> (f a, a)
scanL = undefined -- for now

-- Prefix (left) sums
sumsL :: (Traversable f, Num a) => (a, f a) -> (f a, a)
sumsL = (fmap getSum *** getSum) . scanL . (Sum *** fmap Sum)

counts :: forall f a. (Traversable f, Applicative f, Num a) => (f a, a)
counts = sumsL (0, pure 1 :: f a)

cross :: (Functor g, Functor f) => g a -> f b -> g (f (a , b))
cross as bs = fmap (\ a -> fmap (\ b -> (a,b)) bs) as

dot :: (Applicative f, Foldable f, Num a) => f a -> f a -> a
u `dot` v = sum (liftA2 (*) u v)

-- uroots :: f (Complex t)
-- uroots = fmap (\ k -> exp (-
--  where
--    (indices,n) = counts

ci :: RealFloat t => Complex t
ci = 0 :+ 1

uroot :: RealFloat t => Int -> Complex t
uroot n = exp (- 2 * pi * ci / fromIntegral n)
