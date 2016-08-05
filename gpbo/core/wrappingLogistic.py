##
# wrapping: A program making it easy to use hyperparameter
# optimization software.
# Copyright (C) 2013 Katharina Eggensperger and Matthias Feurer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/usr/bin/env python

from gzip import GzipFile as gfile
import numpy
import time
import sys
import os

import theano
import theano.tensor as T
import collections
import skdata.mnist.dataset

import data_util

__authors__ = ["Jasper Snoek", "Katharina Eggensperger", "Matthias Feurer"]
__contact__ = "automl.org"


class LogisticRegression(object):
    """Multi-class Logistic Regression Class

    The logistic regression is fully described by a weight matrix :math:`W` 
    and bias vector :math:`b`. Classification is done by projecting data 
    points onto a set of hyperplanes, the distance to which is used to 
    determine a class membership probability. 
    """

    def __init__(self, input_, n_in, n_out):
        """ Initialize the parameters of the logistic regression

        :type input_: theano.tensor.TensorType
        :param input_: symbolic variable that describes the input of the
                      architecture (one minibatch)
        
        :type n_in: int
        :param n_in: number of input units, the dimension of the space in 
                     which the datapoints lie

        :type n_out: int
        :param n_out: number of output units, the dimension of the space in 
                      which the labels lie

        """ 

        # initialize with 0 the weights W as a matrix of shape (n_in, n_out) 
        self.W = theano.shared(value=numpy.zeros((n_in, n_out), dtype=
                               theano.config.floatX), name='W')
        # initialize the baises b as a vector of n_out 0s
        self.b = theano.shared(value=numpy.zeros((n_out,), dtype=
                               theano.config.floatX), name='b')


        # compute vector of class-membership probabilities in symbolic form
        self.p_y_given_x = T.nnet.softmax(T.dot(input_, self.W)+self.b)

        # compute prediction as class whose probability is maximal in 
        # symbolic form
        self.y_pred=T.argmax(self.p_y_given_x, axis=1)

        # parameters of the model
        self.params = [self.W, self.b]

        self.L1 = abs(self.W).sum()

        self.L2_sqr = (self.W**2).sum()

    def negative_log_likelihood(self, y):
        """Return the mean of the negative log-likelihood of the prediction
        of this model under a given target distribution.

        .. math::

            \frac{1}{|\mathcal{D}|} \mathcal{L} (\theta=\{W,b\}, \mathcal{D}) =
            \frac{1}{|\mathcal{D}|} \sum_{i=0}^{|\mathcal{D}|} \log(P(Y=y^{(i)}|x^{(i)}, W,b)) \\
                \ell (\theta=\{W,b\}, \mathcal{D})

        :type y: theano.tensor.TensorType
        :param y: corresponds to a vector that gives for each example the
                  correct label

        Note: we use the mean instead of the sum so that
              the learning rate is less dependent on the batch size
        """
        # y.shape[0] is (symbolically) the number of rows in y, i.e., number of examples (call it n) in the minibatch
        # T.arange(y.shape[0]) is a symbolic vector which will contain [0,1,2,... n-1]
        # T.log(self.p_y_given_x) is a matrix of Log-Probabilities (call it LP) with one row per example and one column per class 
        # LP[T.arange(y.shape[0]),y] is a vector v containing [LP[0,y[0]], LP[1,y[1]], LP[2,y[2]], ..., LP[n-1,y[n-1]]]
        # and T.mean(LP[T.arange(y.shape[0]),y]) is the mean (across minibatch examples) of the elements in v,
        # i.e., the mean log-likelihood across the minibatch.
        return -T.mean(T.log(self.p_y_given_x)[T.arange(y.shape[0]),y])

    def errors(self, y):
        """Return a float representing the number of errors in the minibatch 
        over the total number of examples of the minibatch ; zero one
        loss over the size of the minibatch

        :type y: theano.tensor.TensorType
        :param y: corresponds to a vector that gives for each example the 
                  correct label
        """

        # check if y has same dimension of y_pred 
        if y.ndim != self.y_pred.ndim:
            raise TypeError('y should have the same shape as self.y_pred', 
                ('y', target.type, 'y_pred', self.y_pred.type))
        # check if y is of the correct datatype        
        if y.dtype.startswith('int'):
            # the T.neq operator returns a vector of 0s and 1s, where 1
            # represents a mistake in prediction
            return T.mean(T.neq(self.y_pred, y))
        else:
            raise NotImplementedError()

def logistic(params, **kwargs):
    learning_rate = float(params["lrate"])
    l1_reg = 0  # float(params["l1_reg"])
    l2_reg = float(params["l2_reg"])
    batch_size = int(float(params["batchsize"]))
    n_epochs = int(float(params["n_epochs"]))
    
    # This part is taken from logistic_sgd.py
    # but changed to not load data on its own
    
    """
    Demonstrate stochastic gradient descent optimization of a log-linear 
    model

    This is demonstrated on MNIST.
    
    :type learning_rate: float
    :param learning_rate: learning rate used (factor for the stochastic 
                          gradient)

    :type n_epochs: int
    :param n_epochs: maximal number of epochs to run the optimizer 

    :type dataset: string
    :param dataset: the path of the MNIST dataset file from 
                         http://www.iro.umontreal.ca/~lisa/deep/data/mnist/mnist.pkl.gz

    """
    
    if "train" in kwargs and "valid" in kwargs:
        train = kwargs['train']
        valid = kwargs['valid']
        train_targets = kwargs['train_targets']
        valid_targets = kwargs['valid_targets']
    else:
        raise Exception("No training and valid data found")
    
    train_set_x = theano.shared(numpy.asarray(train, dtype=theano.config.floatX))
    train_set_y = theano.shared(numpy.asarray(train_targets, dtype=theano.config.floatX))
    train_set_y = T.cast(train_set_y, 'int32')
    valid_set_x = theano.shared(numpy.asarray(valid, dtype=theano.config.floatX))
    valid_set_y = theano.shared(numpy.asarray(valid_targets, dtype=theano.config.floatX))
    valid_set_y = T.cast(valid_set_y, 'int32')
    print "####", type(valid_set_y), type(train_set_y)
    # compute number of minibatches for training, validation and testing
    n_train_batches = train_set_x.get_value(borrow=True).shape[0] / batch_size
    n_valid_batches = valid_set_x.get_value(borrow=True).shape[0] / batch_size

    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar()    # index to a [mini]batch 
    x     = T.matrix('x')  # the data is presented as rasterized images
    y     = T.ivector('y') # the labels are presented as 1D vector of 
                           # [int] labels

    # construct the logistic regression class
    # Each MNIST image has size 28*28
    classifier = LogisticRegression(input_=x, n_in=28*28, n_out=10)

    # the cost we minimize during training is the negative log likelihood of 
    # the model in symbolic format
    cost = classifier.negative_log_likelihood(y) + l1_reg * classifier.L1 \
         + l2_reg * classifier.L2_sqr
    
    # compiling a Theano function that computes the mistakes that are made by 
    # the model on a minibatch
    """
    test_model = theano.function(inputs = [index], 
            outputs = classifier.errors(y),
            givens={
                x:test_set_x[index*batch_size:(index+1)*batch_size],
                y:test_set_y[index*batch_size:(index+1)*batch_size]})
    """
    validate_model = theano.function(inputs=[index],
            outputs=classifier.errors(y),
            givens={
                x: valid_set_x[index*batch_size:(index+1)*batch_size],
                y: valid_set_y[index*batch_size:(index+1)*batch_size]})

    # compute the gradient of cost with respect to theta = (W,b) 
    g_W = T.grad(cost = cost, wrt = classifier.W)
    g_b = T.grad(cost = cost, wrt = classifier.b)

    # specify how to update the parameters of the model as a dictionary
    updates =collections.OrderedDict(
            [(classifier.W, classifier.W - learning_rate*g_W),
             (classifier.b, classifier.b - learning_rate*g_b)])

    # compiling a Theano function `train_model` that returns the cost, but in 
    # the same time updates the parameter of the model based on the rules 
    # defined in `updates`
    train_model = theano.function(inputs=[index],
                                  outputs=cost,
                                  updates=updates,
                                  givens={
                x:train_set_x[index*batch_size:(index+1)*batch_size],
                y:train_set_y[index*batch_size:(index+1)*batch_size]})

    ###############
    # TRAIN MODEL #
    ###############
    print '... training the model'
    # early-stopping parameters
    patience              = 5000   # look as this many examples regardless
    patience_increase     = 2      # wait this much longer when a new best is
                                   # found
    improvement_threshold = 0.995  # a relative improvement of this much is
                                   # considered significant
    validation_frequency = min(n_train_batches, patience/2)
                                   # go through this many
                                   # minibatche before checking the network
                                   # on the validation set; in this case we
                                   # check every epoch

    best_params          = None
    best_validation_loss = numpy.inf
    start_time = time.clock()

    done_looping = False 
    epoch = 0  
    while (epoch < n_epochs) and (not done_looping):
        epoch += 1
        for minibatch_index in xrange(n_train_batches):

            minibatch_avg_cost = train_model(minibatch_index)
            # iteration number
            iter = epoch * n_train_batches + minibatch_index

            if (iter+1) % validation_frequency == 0: 
                # compute zero-one loss on validation set 
                validation_losses = [validate_model(i) for i in xrange(n_valid_batches)]
                this_validation_loss = numpy.mean(validation_losses)

                print('epoch %i, minibatch %i/%i, validation error %f %%' % \
                    (epoch, minibatch_index+1,n_train_batches, \
                    this_validation_loss*100.))


                # if we got the best validation score until now
                if this_validation_loss < best_validation_loss:
                    #improve patience if loss improvement is good enough
                    if this_validation_loss < best_validation_loss *  \
                       improvement_threshold :
                        patience = max(patience, iter * patience_increase)

                    best_validation_loss = this_validation_loss
                    # test it on the test set

#                    test_losses = [test_model(i) for i in xrange(n_test_batches)]
#                    test_score  = numpy.mean(test_losses)

#                    print(('     epoch %i, minibatch %i/%i, test error of best ' 
#                       'model %f %%') % \
#                        (epoch, minibatch_index+1, n_train_batches,test_score*100.))
                    print ('     epoch %i, minibatch %i/%i') % \
                          (epoch, minibatch_index+1, n_train_batches)
            #if patience <= iter :
            #    done_looping = True
            #    break

    end_time = time.clock()
#    print(('Optimization complete with best validation score of %f %%,'
#           'with test performance %f %%') %  
#                 (best_validation_loss * 100., test_score*100.))
    print(('Optimization complete with best validation score of %f %%') %  
                 (best_validation_loss * 100.,))    

    print 'The code run for %d epochs, with %f epochs/sec' % \
          (epoch, 1. * epoch / (end_time - start_time))
    print >> sys.stderr, ('The code for file ' + os.path.split(__file__)[1] +
                          ' ran for %.1fs' % ((end_time - start_time)))
    return best_validation_loss


def main(params_, **kwargs):
    print 'Params: ', params_
    #raise IndexError('aaaaaaaaaaaaaaaaaaaaaaaaaaaaa!')
    dataset = skdata.mnist.dataset.MNIST()
    dataset.fetch(True)
    convex_inputs = skdata.mnist.dataset.read(gfile(os.path.join(dataset.home(),"train-images-idx3-ubyte.gz")))
    convex_labels = skdata.mnist.dataset.read(gfile(os.path.join(dataset.home(), "train-labels-idx1-ubyte.gz")))
    convex_test = skdata.mnist.dataset.read(gfile(os.path.join(dataset.home(), "t10k-images-idx3-ubyte.gz")))
    convex_test_labels = skdata.mnist.dataset.read(gfile(os.path.join(dataset.home(), "t10k-labels-idx1-ubyte.gz")))
    convex_inputs = numpy.array(convex_inputs.reshape((-1, 784)), dtype=numpy.float32)
    convex_test = numpy.array(convex_test.reshape((-1, 784)), dtype=numpy.float32)

    fold = kwargs['fold']
    folds = kwargs['folds']
    downsize=kwargs['downsize']
    if folds == 1:
        train = convex_inputs[:int(50000*downsize)]
        valid = convex_inputs[50000:60000]
        test = convex_test
        train_targets = convex_labels[:int(50000*downsize)]
        valid_targets = convex_labels[50000:60000]
        test_targets = convex_test_labels
    elif folds > 1:
        cv_data = numpy.copy(convex_inputs[60000])
        train, valid = data_util.prepare_cv_for_fold(cv_data, fold, folds)
        cv_labels = numpy.copy(convex_labels[:60000])
        train_targets, valid_targets = data_util.prepare_cv_for_fold(cv_labels,
                                                                     fold, folds)
        test = convex_test
        test_targets = convex_test_labels
    else:
        raise ValueError("Folds cannot be less than 1")

    # I do not know why 256, but this gives us equal data than the one obtained
    # by Jasper Snoek
    train /= 256
    valid /= 256

    kwargs['train'] = train
    kwargs['train_targets'] = train_targets
    kwargs['valid'] = valid
    kwargs['valid_targets'] = valid_targets
    y = logistic(params_, **kwargs)
    print 'Result: ', y
    return y
