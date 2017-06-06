import tensorflow as tf
import numpy as np
import GPflow as gpf
from GPflow.param import Param
from GPflow import transforms

DXmax = 2

np_float_type = gpf.kernels.np_float_type


class testk(gpf.kernels.Stationary):
    def K(self, X, Y=None, presliced=False):
        if not presliced:
            X, Y = self._slice(X, Y)
        if Y is None:
            Y = X

        X = X / self.lengthscales
        Y = Y / self.lengthscales
        K = self.variance * tf.exp(-0.5 * tf.square(X - tf.transpose(Y)))
        K_ = tf.Print(K, [K], summarize=100, message='\n\nk')
        return K_


class Stationary(gpf.kernels.Kern):
    """
    Stationary kernel with derivatives incicated my auxiliary dimensions
    """

    def __init__(self, input_dim, **kwargs):
        gpf.kernels.Kern.__init__(self, 2 * input_dim, **kwargs)
        self.latent_input_dim = input_dim
        if not 'lengthscales' in kwargs.keys():
            lengthscales = np.ones(input_dim, np_float_type)
        else:
            # accepts float or array:
            lengthscales = kwargs['lengthscales'] * np.ones(input_dim, np_float_type)
        self.lengthscales = Param(lengthscales, transforms.positive)

        if not 'variance' in kwargs.keys():
            variance = 1.
        self.variance = Param(variance, transforms.positive)
    def Rsplit(self, X, Y):
        """
        form a 3tensor of signed distances [nx * ny * d] per dimension and matrix of squaredeuclid distances [nx * ny]
        """
        D = tf.expand_dims(X, axis=1) - tf.expand_dims(Y, axis=0)
        Rs = tf.reduce_sum(tf.square(D), axis=2)
        # D = tf.Print(D,[tf.shape(D)],summarize=100,message="yyy")
        return D, Rs

    def polyconv(self, M):
        """
        compute the coefficients of the product of polynomials with individual coefs as the rows of M. Fixed width
        based on DXmax
        """
        w = 4 * DXmax + 1
        N = tf.concat([tf.zeros([self.latent_input_dim, 4]), tf.cast(M, tf.float32)], axis=1)
        P = tf.foldl(lambda a, x: tf.reshape(
            tf.nn.convolution(tf.reshape(a, [1, w, 1]), tf.reshape(tf.reverse_v2(x, axis=[0]), [w, 1, 1]), 'SAME'),
            [w]), N)
        Q = tf.slice(P, [4], [-1])
        # Q=tf.Print(Q,[M,Q,tf.shape(Q)],summarize=100,message='q')
        return tf.cast(Q, tf.float64)

    def condense(self, M, D):
        I = tf.ones_like(D)

        Eq = [tf.cast(tf.equal(M, i), tf.float64) for i in range(5)]
        H = tf.stack([Eq[0] * I,
                      Eq[1] * D * 2 + Eq[2] * 2 * I,
                      Eq[2] * 4 * D ** 2 + Eq[3] * 12 * D + Eq[4] * 12 * I,
                      Eq[3] * 8 * D ** 3 + Eq[4] * 48 * D ** 2,
                      Eq[4] * 16 * D ** 4], axis=3)
        pad = tf.slice(tf.zeros_like(H), [0, 0, 0, 1], [-1, -1, -1, -1])
        A = tf.map_fn(lambda x: tf.map_fn(self.polyconv, x), H)
        return A

    def Dmap(self, Dx, Dy):
        """
        form a 3tensor of derivative indicies [nx * ny * d] per dimension and matrix of total derivative [nx * ny]
        """
        M = tf.expand_dims(tf.cast(Dx, tf.int32), axis=1) + tf.expand_dims(tf.cast(Dy, tf.int32), axis=0)
        T = tf.reduce_sum(M, axis=2)
        return M, T

    def dZdX(self, Dx, Dy):
        DyFlat = tf.reduce_sum(Dy, axis=1, keep_dims=True)
        DxFlat = tf.reduce_sum(Dx, axis=1, keep_dims=True)
        DRX = tf.reduce_prod(tf.pow(tf.cast(self.lengthscales, tf.float64), -Dx), axis=1,
                             keep_dims=True) * tf.transpose(tf.ones_like(DyFlat))
        DRY = tf.ones_like(DxFlat) * tf.transpose(
            tf.reduce_prod(tf.pow(-tf.cast(self.lengthscales, tf.float64), -Dy), axis=1, keep_dims=True))
        DR = tf.multiply(DRX, DRY)
        return DR

    def Kdiag(self, X):
        return tf.squeeze(tf.diag_part(self.K(X)))


class dkern(Stationary):
    def K(self, X_, Y_=None, presliced=False):
        if not presliced:
            X_, Y_ = self._slice(X_, Y_)

        if Y_ is None:
            Y_ = X_

        X, Dx = tf.split(X_, num_or_size_splits=2, axis=1)
        Y, Dy = tf.split(Y_, num_or_size_splits=2, axis=1)
        Zx = X / self.lengthscales
        Zy = Y / self.lengthscales

        M, T = self.Dmap(Dx, Dy)
        D, Rs = self.Rsplit(Zx, Zy)

        K = self.variance * tf.exp(-0.5 * Rs)
        dKdU = tf.stack([K, -0.5 * K, 0.25 * K, -0.125 * K, 0.0625 * K], axis=2)

        dUdZ = self.condense(M, D)
        dZdX = self.dZdX(Dx, Dy)

        Kpost = tf.reduce_sum(dKdU * dUdZ, axis=2) * dZdX
        return Kpost

class Pointwise_Hetroskedastic(gpf.kernels.Kern):
    """
    The White kernel
    """
    def __init__(self, input_dim, active_dims=None):
        gpf.kernels.Kern.__init__(self, input_dim, active_dims)

    def K(self, X, X2=None, presliced=False):
        if not presliced:
            X, X2 = self._slice(X,X2)
        if X2 is None:
#            X = tf.Print(X,[X],message='x',summarize=100)
            d = tf.zeros(tf.stack([tf.shape(X)[0]]),gpf.kernels.float_type)
            return tf.diag(tf.squeeze(X))
        else:
            shape = tf.stack([tf.shape(X)[0], tf.shape(X2)[0]])
            return tf.zeros(shape,gpf.kernels.float_type)
    def Kdiag(self,X):
        K= tf.squeeze(tf.diag_part(self.K(X)))
        return K

