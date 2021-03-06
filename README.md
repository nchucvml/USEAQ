# USEAQ
## USEAQ: Ultra-fast Superpixel Extraction via Adaptive Sampling from Quantized Regions

Abstract — We present a novel and highly efficient superpixel extraction method called USEAQ to generate regular and compact
superpixels in an image. To reduce the computational cost of iterative optimization procedures adopted in most recent
approaches, the proposed USEAQ for superpixel generation works in a one-pass fashion. It firstly performs joint spatial
and color quantizations and groups pixels into regions. It then takes into account the variations between regions, and adaptively
samples one or a few superpixel candidates for each region. It finally employs maximum a posteriori (MAP) estimation to
assign pixels to the most spatially consistent and perceptually similar superpixels. It turns out that the proposed USEAQ is
quite efficient, and the extracted superpixels can precisely adhere to boundaries of objects. Experimental results show that USEAQ
achieves better or equivalent performance compared to the state-of-the-art superpixel extraction approaches in terms of boundary
recall, undersegmentation error, achievable segmentation accuracy, the average miss rate, average undersegmentation error,
and average unexplained variation, and it is significantly faster than these approaches.

* Flowchart of USEAQ

![image](https://github.com/nchucvml/USEAQ/blob/master/Fig_Flowchart.png)

* If you use the USEAQ or USEQ codes (USEQsuperpixel_ICPR.zip) for your research, please cite the following papers.

[1] Chun-Rong Huang, Wei-Cheng Wang, Wei-An Wang, Szu-Yu Lin, and Yen-Yu Lin, "USEAQ: Ultra-fast Superpixel Extraction via Adaptive Sampling from Quantized Regions,"  IEEE Transactions on Image Processing, vol. 27, no. 10, pp. 4916-4931, Oct. 2018. 

[2] Chun-Rong Huang, Wei-An Wang, Szu-Yu Lin, and Yen-Yu Lin, "USEQ: Ultra-Fast Superpixel Extraction via Quantization," in Proc. International Conference on Pattern Recognition, ICPR’16, pp. 1966-1971, Dec. 2016.

* The webpage of USEQ is at http://cvml.cs.nchu.edu.tw/USEQ.html 

As for the evaluation of Stanford Background Dataset, please refer to https://davidstutz.de/projects/superpixel-benchmark/

As for the evaluation of BSDS500 dataset, please refer to https://www.tu-chemnitz.de/etit/proaut/en/research/superpixel.html

