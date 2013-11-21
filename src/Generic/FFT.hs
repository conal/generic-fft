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

import Data.Functor ((<$>))
import Data.Foldable (sum)
import Data.Traversable (Traversable(..), mapAccumL)
import Control.Applicative (Applicative(..))
import Data.Monoid (Monoid(..),Product(..))

import Data.Complex (Complex(..))

import TypeUnary.Nat

import Data.UniformPair

import qualified Data.FTree.BottomUp as I
import qualified Data.FTree.TopDown  as O

-- | Input trees -- bottom-up
type I = I.T Pair

-- | Output trees -- bottom-up
type O = O.T Pair

-- | FFT computation, from one functor to another
class HasFFT f g | f -> g where
  fft :: RealFloat t => f (Complex t) -> g (Complex t)

instance HasFFT Pair Pair where
  fft (a :# b) = a+b :# a-b

instance IsNat n => HasFFT (I n) (O n) where
  fft = fft' nat

{-
-- fft' :: RealFloat t => Nat n -> I n (Complex t) -> O n (Complex t)
fft' :: (Traversable f, Applicative f, HasFFT f f, RealFloat t) =>
        Nat n -> (f :+^ n) (Complex t) -> (f :^+ n) (Complex t)
fft' Zero     = O.L . I.unL
fft' (Succ m) =
  O.B . (inTranspose.fmap) fft . tweakParts . fmap (fft' m) . transpose . I.unB

-- tweakParts :: Unop (Pair (O n (Complex t)))
-- tweakParts = secondP undefined -- for now
tweakParts :: Unop (f ((f :^+ n) (Complex t)))
tweakParts = undefined -- for now

-}

-- type (:+^) = I.T    -- bottom-up
-- type (:^+) = O.T    -- top-down

-- fft' :: (Traversable f, Traversable g, Applicative f, Applicative g, HasFFT f g, RealFloat t) =>
--         Nat n -> (f :+^ n) (Complex t) -> (g :^+ n) (Complex t)

fft' :: (Traversable f, Traversable g, Applicative f, Applicative g, HasFFT f g, RealFloat t) =>
        Nat n -> I.T f n (Complex t) -> O.T g n (Complex t)

fft' Zero     = O.L . I.unL
fft' (Succ m) =
  O.B . (inTranspose.fmap) fft . {- tweakParts . -} fmap (fft' m) . transpose . I.unB


fftC :: ( Traversable f, Applicative f, Applicative g, Traversable h, Traversable k, Applicative k
        , HasFFT k g, HasFFT h f
        , RealFloat t ) =>
        h (k (Complex t)) -> g (f (Complex t))
fftC = transpose . fmap fft . transpose . id . fmap fft . transpose

-- fftC = (inTranspose.fmap) fft . id . fmap fft . transpose


-- fftC :: -- (Traversable f, Traversable g, Applicative f, Applicative g, HasFFT f g, RealFloat t) =>
--         g (f (Complex t)) -> (g :^+ n) (Complex t)
-- fftC = O.L . I.unL
-- fftC =
--   O.B . (inTranspose.fmap) fft . {- tweakParts . -} fmap (fft' m) . transpose . I.unB

-- I could also drop the initial and/or final transposition, e.g.,

-- fftC' :: ( Traversable f, Applicative f, Applicative g, Traversable h, Traversable k, Applicative k
--         , HasFFT k g, HasFFT h f
--         , RealFloat t ) =>
--         k (h (Complex t)) -> f (g (Complex t))

fftC' :: (Applicative f', Traversable g, HasFFT g g', HasFFT f f', RealFloat t) =>
         g (f (Complex t)) -> f' (g' (Complex t))
fftC' = fmap fft . transpose . twiddle . fmap fft

-- fmap fft  :: g  (f  c) -> g  (f' c)
-- twiddle   :: g  (f' c) -> g  (f' c)
-- transpose :: g  (f' c) -> f' (g  c)
-- fmap fft  :: f' (g  c) -> f' (g' c)

fftC'' :: ( Traversable f, Traversable g, Applicative f, Applicative g'
          , HasFFT f f', HasFFT g g', RealFloat t ) =>
          g (f (Complex t)) -> g' (f' (Complex t))
fftC'' = fmap fft . transpose . twiddle . fmap fft . transpose

-- transpose :: g  (f  c) -> f  (g  c)
-- fmap fft  :: f  (g  c) -> f  (g' c)
-- twiddle   :: f  (g' c) -> f  (g' c)
-- transpose :: f  (g' c) -> g' (f  c)
-- fmap fft  :: g' (f  c) -> g' (f' c)

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

scanL :: (Traversable f, Monoid a) => (a, f a) -> (f a, a)
scanL = undefined -- for now

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
