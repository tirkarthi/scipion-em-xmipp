import pyworkflow.em.metadata as md
import xmippLib
import numpy as np
import random

from sklearn.utils import shuffle
from sklearn.cross_validation import train_test_split
from .augmentators import (_random_flip_leftright, _random_flip_updown, _mismatch_projection,
                          _random_90degrees_rotation, _random_rotation,generateEmptyParticlesFunction)


BATCH_SIZE= 16
def getDataGenerator( imgsMdXmd, masksMdXmd, xmdEmptyParts=None, augmentData=True, nEpochs=-1, isTrain=True, valFraction=0.1, 
                      batchSize= BATCH_SIZE, doTanhNormalize=False, simulateEmptyParts=True, addMismatch=False): 

  if nEpochs<1: 
    nEpochs= 9999999
  mdImgs  = md.MetaData(imgsMdXmd)
  mdMasks = md.MetaData(masksMdXmd)
    
  nImages= int(mdImgs.size())

  I= xmippLib.Image()      

  imgFnames = mdImgs.getColumnValues(md.MDL_IMAGE)
  maskFnames= mdMasks.getColumnValues(md.MDL_IMAGE)
  
  I.read( imgFnames[0] )
  shape= I.getData().shape+ (1,)

  if not xmdEmptyParts is None:
    mdEmpty= md.MetaData(xmdEmptyParts)
    nImages+= int(mdEmpty.size())
    emptyFnames= mdEmpty.getColumnValues(md.MDL_IMAGE)
    imgFnames+= emptyFnames
    maskFnames+= [None]*len(emptyFnames)
  
  stepsPerEpoch= nImages//batchSize
    

  if augmentData:
    augmentFuns= [_random_flip_leftright, _random_flip_updown, _random_90degrees_rotation, _random_rotation ]
    if simulateEmptyParts==True:
      augmentFuns+= [generateEmptyParticlesFunction(shape, prob=0.2)]
    if addMismatch==True:
      augmentFuns+= [_mismatch_projection]
    def augmentBatch( batchX, batchY):
      for fun in augmentFuns:
        if bool(random.getrandbits(1)):
          batchX, batchY= fun(batchX, batchY)
      return batchX, batchY
  else:
    def augmentBatch( batchX, batchY): return batchX, batchY
    
  if valFraction>0:
    (imgFnames_train, imgFnames_val, maskFnames_train,
     maskFnames_val) = train_test_split(imgFnames, maskFnames, test_size=valFraction, random_state=121)  
    if isTrain:
      imgFnames, maskFnames= imgFnames_train, maskFnames_train
    else:
      imgFnames, maskFnames= imgFnames_val, maskFnames_val
 
  if doTanhNormalize:
    def readImgAndMask(fnImageImg, fnImageMask):
      I.read(fnImageImg)
      img= normalization( np.expand_dims(I.getData(), -1), use0_1_norm_instead_1_1=False)
      if fnImageMask is None:
        mask= -1*np.ones_like(img)
      else:
        I.read(fnImageMask)
        mask= normalization(np.expand_dims(I.getData(), -1), use0_1_norm_instead_1_1=False)
      return img, mask
  else:
    def readImgAndMask(fnImageImg, fnImageMask):
      I.read(fnImageImg)
      img=np.expand_dims(I.getData(),-1)
      if fnImageMask is None:
        mask= -1*np.ones_like(img)
      else:
        I.read(fnImageMask)
      mask= np.expand_dims(I.getData(),-1)
      return img, mask
      
  def dataIterator(imgFnames, maskFnames, nBatches=None):
    
    batchStack = np.zeros((batchSize,)+shape )
    batchLabels = np.zeros((batchSize,)+shape )
    currNBatches=0 
    for epoch in range(nEpochs):
      if isTrain:
        imgFnames, maskFnames= shuffle(imgFnames, maskFnames)
      n=0
      for fnImageImg, fnImageMask in zip(imgFnames, maskFnames):
        batchStack[n,...], batchLabels[n,...]= readImgAndMask(fnImageImg, fnImageMask)
        n+=1
        if n>=batchSize:
          yield augmentBatch(batchStack, batchLabels)
          n=0
          currNBatches+=1
          if nBatches and currNBatches>=nBatches:
            break
      if n>0:
        yield augmentBatch(batchStack[:n,...], batchLabels[:n,...])
        
  return dataIterator(imgFnames, maskFnames), stepsPerEpoch
  
def extractNBatches(valIterator, maxBatches=-1):
  x_val=[]
  y_val=[]
  for i, (x, y) in enumerate(valIterator):
    x_val.append(x)
    y_val.append(y)
    if i== maxBatches:
      break
  return ( np.concatenate(x_val, axis=0), np.concatenate(y_val, axis=0 ))

def normalization( img, use0_1_norm_instead_1_1=True):
  normData= (img -np.min(img))/ (np.max(img)-np.min(img))
  if not use0_1_norm_instead_1_1:
    normData= 2*normData -1
  if np.any( np.isnan(normData)):
    normData= np.zeros_like(normData)
  return normData
  
def normalizeImgs(batch_img, use0_1_norm_instead_1_1=True):
  for i in range(batch_img.shape[0]):
    batch_img[i]= normalization(batch_img[i], use0_1_norm_instead_1_1)
  return batch_img