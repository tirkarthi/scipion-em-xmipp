import sys
import keras
import math
from .DeepLearningGeneric import DeepLearningModel
from .dataGenerator import getDataGenerator, extractNBatches, BATCH_SIZE

NUM_BATCHES_PER_EPOCH= 50

class UNET(DeepLearningModel):

  
  def __init__(self, boxSize, saveModelDir, gpuList="0", batchSize=BATCH_SIZE, init_n_filters=64):
  
    DeepLearningModel.__init__(self,boxSize, saveModelDir, gpuList, batchSize)
    self.epochSize= NUM_BATCHES_PER_EPOCH
    self.init_n_filters= init_n_filters   
           
  def _setShape(self, boxSize):
    self.shape= (boxSize,boxSize,1)
    return self.shape
    
  def train(self, learningRate, nEpochs, xmdParticles, xmdProjections):        
    try:
      print("loading model ")
      model = keras.models.load_model(self.saveModelDir, custom_objects={})
      print("previous model loaded")
    except Exception as e:
      print(e)
      model = self.UNet( img_shape=self.shape, start_ch=self.init_n_filters, batchnorm=False )
      optimizer= keras.optimizers.Adam(lr= learningRate, beta_1=0.9, beta_2=0.999, epsilon=1e-8, decay=0)
      model.compile(loss='mse', optimizer=optimizer)
      
      reduceLrCback= keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=3, verbose=1, mode='auto', min_lr=1e-8)
      earlyStopCBack= keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0, patience=5, verbose=1, mode='auto')
      saveModel= keras.callbacks.ModelCheckpoint(self.saveModelDir, monitor='val_loss', verbose=1, save_best_only=True)

      trainIterator, stepsPerEpoch= getDataGenerator(xmdParticles, xmdProjections, isTrain=True, valFraction=0.1)
      
      print("train/val split"); sys.stdout.flush()
      valIterator, stepsPerEpoch_val= getDataGenerator(xmdParticles, xmdProjections, augmentData=False, 
                                                             isTrain=False, valFraction=0.1)
      valData= extractNBatches(valIterator, min(10, stepsPerEpoch_val))
      del valIterator
#      valData=None
      nEpochs_init= nEpochs
      nEpochs= max(1, nEpochs_init*float(stepsPerEpoch)/self.epochSize)
      print("nEpochs : %.1f --> Epochs: %d.\nTraining begins: Epoch 0/%d"%(nEpochs_init, nEpochs, nEpochs))
      sys.stdout.flush()
      
      model.fit_generator(trainIterator, epochs= nEpochs, steps_per_epoch=self.epochSize, #steps_per_epoch=stepsPerEpoch,
                        verbose=2, callbacks=[reduceLrCback, earlyStopCBack, saveModel],
                        validation_data=valData, max_queue_size=12, workers=1, use_multiprocessing=False)
                        
      
  def conv_block(self, m, dim, acti, bn, res, do=0):
      n = keras.layers.Conv2D(dim, 3, activation=acti, padding='same')(m)
      n = keras.layers.BatchNormalization()(n) if bn else n
      n = keras.layers.Dropout(do)(n) if do else n
      n = keras.layers.Conv2D(dim, 3, activation=acti, padding='same')(n)
      n = keras.layers.BatchNormalization()(n) if bn else n
      return keras.layers.Concatenate()([m, n]) if res else n

  def level_block(self, m, dim, depth, inc, acti, do, bn, mp, up, res):
      if depth > 0:
          n = self.conv_block(m, dim, acti, bn, res)
          m = keras.layers.MaxPooling2D()(n) if mp else keras.layers.Conv2D(dim, 3, strides=2, padding='same')(n)
          m = self.level_block(m, int(inc*dim), depth-1, inc, acti, do, bn, mp, up, res)
          if up:
              m = keras.layers.UpSampling2D()(m)
              m = keras.layers.Conv2D(dim, 2, activation=acti, padding='same')(m)
          else:
              m = keras.layers.Conv2DTranspose(dim, 3, strides=2, activation=acti, padding='same')(m)
          n = keras.layers.Concatenate()([n, m])
          m = self.conv_block(n, dim, acti, bn, res)
      else:
          m = self.conv_block(m, dim, acti, bn, res, do)
      return m

  def UNet(self, img_shape, out_ch=1, start_ch=32, depth=3, inc_rate=2., activation='relu', 
      dropout=0.5, batchnorm=False, maxpool=True, upconv=True, residual=False):
      
      i = keras.layers.Input(shape=img_shape)
      PAD_SIZE= ( 2**(int(math.log(img_shape[1], 2))+1)-img_shape[1]) //2
      x = keras.layers.ZeroPadding2D( (PAD_SIZE, PAD_SIZE) )( i ) #Padded to 2**N
      o = self.level_block(x, start_ch, depth, inc_rate, activation, dropout, batchnorm, maxpool, upconv, residual)
      o = keras.layers.Conv2D(out_ch, 1, activation='linear')(o)
      o = keras.layers.Lambda(lambda m: m[:,PAD_SIZE:-PAD_SIZE,PAD_SIZE:-PAD_SIZE,:] )( o )
      return  keras.models.Model(inputs=i, outputs=o)


def segmentationLoss(trn_labels_batch, logits):
  '''
  '''
  logits=tf.reshape(logits, [-1])
  trn_labels=tf.reshape(trn_labels_batch, [-1])
  inter=tf.reduce_sum(tf.multiply(logits,trn_labels))
  union=tf.reduce_sum(tf.subtract(tf.add(logits,trn_labels),tf.multiply(logits,trn_labels)))
  loss=tf.subtract(tf.constant(1.0, dtype=tf.float32),tf.div(inter,union))
  return loss



#####################################################################################################################
# AltModelCheckpoint taken from https://github.com/TextpertAi/alt-model-checkpoint/blob/master/alt_model_checkpoint/__init__.py#L9     #
#####################################################################################################################
from keras.callbacks import ModelCheckpoint


class AltModelCheckpoint(ModelCheckpoint):
    def __init__(self, filepath, alternate_model, **kwargs):
        """
        Additional keyword args are passed to ModelCheckpoint; see those docs for information on what args are accepted.
        :param filepath:
        :param alternate_model: Keras model to save instead of the default. This is used especially when training multi-
                                gpu models built with Keras multi_gpu_model(). In that case, you would pass the original
                                "template model" to be saved each checkpoint.
        :param kwargs:          Passed to ModelCheckpoint.
        """

        self.alternate_model = alternate_model
        super(type(self),self).__init__(filepath, **kwargs)

    def on_epoch_end(self, epoch, logs=None):
        model_before = self.model
        self.model = self.alternate_model
        super(type(self),self).on_epoch_end(epoch, logs)
        self.model = model_before
        
