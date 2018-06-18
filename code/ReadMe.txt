Please contact wang.wei.cheng.tj{at}gmail(dot)com if you encounter any problem

******************************************
Parameter Setting
******************************************
m_fTheta - default as 4
  m_fTheta is the number of each color channel be divided into.

m_fNumCandidates - default as 10
  m_fNumCandidates is maximum number of the adaptive sampling of superpixel candidates S in each spatially quantized region.

para_fOmega - default as 0.01
  Please refer to Section IV-B. Parameter Selection for more detail
  
vint_superpixelNum
  the number of superpixel you wish to generate
  
bParallelBoosting
  Whether to use parallel boosting (in Quantizations and MAP)


******************************************
Output Option
******************************************
op_bStoreLabelImage
  whether to store the superpixel extraction result stored in int16, will generate 1-channel Mat for evaluation

op_bStoreColourImage
  whether to store the superpixel extraction result presented in different colour

op_bStoreBoundedImage
  whether to store the superpixel extraction result presented by drawing the contour on the original image

op_bShowResults
  whether to show results in different colour and contoured image on window

path_strImageFolderPath
  path of image folder, will scan {path of image folder}/*.jpg or {path of image folder}/*.png

path_strResultFolder
  path of result folder
