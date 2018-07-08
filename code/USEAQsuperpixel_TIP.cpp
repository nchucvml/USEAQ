#include "stdafx.h"
#include "USEAQsuperpixel_TIP.h"
#include <vector>
#include <algorithm>
#include <functional>
#define NORMALIZE 255.0f

std::string lazyResultPath;
std::string  resultPath;


USEAQsuperpixel_TIP::USEAQsuperpixel_TIP(void)
	: m_xNum(1), m_yNum(1), m_matData(nullptr), m_imageGrid(nullptr), dominateColorNum(3)
{
	m_fTheta = 4;
	m_fNumCandidates = 0.1f;
	m_fRefinementMag = 1.0f / 2.5f;
}

USEAQsuperpixel_TIP::~USEAQsuperpixel_TIP(void)
{
	if (m_imageGrid != nullptr) { delete[]m_imageGrid; m_imageGrid = nullptr; }
}

bool USEAQsuperpixel_TIP::InitializeSp(cv::Mat &matData, int xNum, int yNum)
{
	m_inputTemp = matData.clone();

	m_matData = &m_inputTemp;

	resultPath = lazyResultPath;

	if (!m_labelData.empty())
		m_labelData.release();
	m_labelData = cv::Mat::zeros(m_matData->size(), CV_32S);
	if (m_imageGrid != nullptr)
		delete[]m_imageGrid; m_imageGrid = nullptr;

	m_xNum = ((xNum > 0 && xNum <= m_matData->cols) ? xNum : m_xNum);
	m_yNum = ((yNum > 0 && yNum <= m_matData->rows) ? yNum : m_yNum);
	m_imageGrid = new CGrid[m_yNum * m_xNum];

	return true;
}

//  BuildImageGridbyQuantization  //
bool USEAQsuperpixel_TIP::BuildImageGridbyQuantization(int begin, int end)
{
	if (m_imageGrid == nullptr || m_matData->empty())
		return false;

	int xPixelNum = m_matData->cols / m_xNum;
	int yPixelNum = m_matData->rows / m_yNum;

	int xRemainNum = m_matData->cols % m_xNum;
	int yRemainNum = m_matData->rows % m_yNum;

	m_xPixelNum = (xRemainNum == 0 ? xPixelNum : xPixelNum + 1) << 1;
	m_yPixelNum = (yRemainNum == 0 ? yPixelNum : yPixelNum + 1) << 1;

	int qtzlv = m_fTheta;
	int colorQtzize = ceil(255.0f / qtzlv);

	begin = (begin == -1 ? 0 : begin);
	end = (end == -1 ? m_yNum : end);

	int qtzTable[8][8][8];
	float qtzColor[8][8][8][3];
	float qtzSpatial[8][8][8][2];

	int dominQtzLabelCount[10];
	for (int i = begin; i < end; i++)
	{
		CGrid *gridPtr = &m_imageGrid[i * m_xNum];

		// calculate quantization labels for each grid
		for (int j = 0; j < m_xNum; j++)
		{
			gridPtr->x = xPixelNum * j + (j < xRemainNum ? j : xRemainNum);
			gridPtr->y = yPixelNum * i + (i < yRemainNum ? i : yRemainNum);

			gridPtr->width = (j < xRemainNum ? xPixelNum + 1 : xPixelNum);
			gridPtr->height = (i < yRemainNum ? yPixelNum + 1 : yPixelNum);

			gridPtr->size = gridPtr->width * gridPtr->height;

			float *meanVectorPtr = gridPtr->meanVector;
			uchar *pixelPtr = (uchar*)m_matData->data + ((gridPtr->y * m_matData->cols + gridPtr->x) * m_matData->channels());

			memset(qtzTable, 0, sizeof(qtzTable));
			memset(qtzColor, 0, sizeof(qtzColor));
			memset(qtzSpatial, 0, sizeof(qtzSpatial));
			memset(dominQtzLabelCount, 0, sizeof(dominQtzLabelCount));


			int r, g, b;
			float colorR, colorG, colorB;


			float n_tempR, n_tempG, n_tempB;
			float n_tempX, n_tempY;



			for (int k = gridPtr->y; k < gridPtr->y + gridPtr->height; k++)
			{
				for (int l = gridPtr->x; l < gridPtr->x + gridPtr->width; l++)
				{
					meanVectorPtr[0] += ((float)l / (float)m_xPixelNum);
					meanVectorPtr[1] += ((float)k / (float)m_yPixelNum);

					r = (int)(*pixelPtr) / colorQtzize;
					colorR = (float)*pixelPtr++;
					g = (int)(*pixelPtr) / colorQtzize;
					colorG = (float)*pixelPtr++;
					b = (int)(*pixelPtr) / colorQtzize;
					colorB = (float)*pixelPtr++;

					qtzTable[r][g][b]++;
					qtzColor[r][g][b][0] = qtzColor[r][g][b][0] + colorR;
					qtzColor[r][g][b][1] = qtzColor[r][g][b][1] + colorG;
					qtzColor[r][g][b][2] = qtzColor[r][g][b][2] + colorB;
					qtzSpatial[r][g][b][0] = qtzSpatial[r][g][b][0] + l;
					qtzSpatial[r][g][b][1] = qtzSpatial[r][g][b][1] + k;
				}


				pixelPtr += (m_matData->cols - gridPtr->width) * m_matData->channels();
			}
			const int SUPSZ = m_matData->rows*m_matData->cols / (m_xNum * m_yNum);
			int maxSize = 0;
			float totalR = 0, totalG = 0, totalB = 0, totalSize = 0;
			int c = 0;

			for (int idxR = 0; idxR < 8; idxR++)
				for (int idxG = 0; idxG < 8; idxG++)
					for (int idxB = 0; idxB < 8; idxB++)
					{
						totalR = qtzColor[idxR][idxG][idxB][0] + totalR;
						totalG = qtzColor[idxR][idxG][idxB][1] + totalG;
						totalB = qtzColor[idxR][idxG][idxB][2] + totalB;

						if (qtzTable[idxR][idxG][idxB] > SUPSZ *m_fNumCandidates) {
							n_tempR = qtzColor[idxR][idxG][idxB][0] / qtzTable[idxR][idxG][idxB];
							n_tempG = qtzColor[idxR][idxG][idxB][1] / qtzTable[idxR][idxG][idxB];
							n_tempB = qtzColor[idxR][idxG][idxB][2] / qtzTable[idxR][idxG][idxB];
							n_tempX = qtzSpatial[idxR][idxG][idxB][0] / qtzTable[idxR][idxG][idxB];
							n_tempY = qtzSpatial[idxR][idxG][idxB][1] / qtzTable[idxR][idxG][idxB];

							gridPtr->dominateColorVector[c][0] = n_tempR / NORMALIZE;
							gridPtr->dominateColorVector[c][1] = n_tempG / NORMALIZE;
							gridPtr->dominateColorVector[c][2] = n_tempB / NORMALIZE;

							gridPtr->dominateSpatialVector[c][0] = n_tempX / (float)m_xPixelNum;
							gridPtr->dominateSpatialVector[c][1] = n_tempY / (float)m_yPixelNum;

							gridPtr->dominateLabelSize[c] = qtzTable[idxR][idxG][idxB];
							gridPtr->dominateLabelCount++;
							c++;
						}

					}
			printf("");

			if (c == 0)
			{
				gridPtr->dominateColorVector[c][0] = (totalR / gridPtr->size) / NORMALIZE;
				gridPtr->dominateColorVector[c][1] = (totalG / gridPtr->size) / NORMALIZE;
				gridPtr->dominateColorVector[c][2] = (totalB / gridPtr->size) / NORMALIZE;
				gridPtr->dominateSpatialVector[c][0] = (gridPtr->x + gridPtr->width / 2.0f) / (float)m_xPixelNum;
				gridPtr->dominateSpatialVector[c][1] = (gridPtr->y + gridPtr->height / 2.0f) / (float)m_yPixelNum;
			}

			gridPtr->meanVector[0] /= (float)(gridPtr->size);
			gridPtr->meanVector[1] /= (float)(gridPtr->size);
			gridPtr->label = i * m_xNum + j;
			gridPtr->updateMeanVector[0] = gridPtr->meanVector[0];
			gridPtr->updateMeanVector[1] = gridPtr->meanVector[1];
			// the normalize color
			gridPtr->updateMeanVector[2] = gridPtr->meanVector[2];
			gridPtr->updateMeanVector[3] = gridPtr->meanVector[3];
			gridPtr->updateMeanVector[4] = gridPtr->meanVector[4];

			gridPtr++;
		}
	}
	return true;
}

bool USEAQsuperpixel_TIP::BuildImageGrid(int begin, int end)
{
	if (m_imageGrid == nullptr || m_matData->empty())
		return false;

	int xPixelNum = m_matData->cols / m_xNum;
	int yPixelNum = m_matData->rows / m_yNum;

	int xRemainNum = m_matData->cols % m_xNum;
	int yRemainNum = m_matData->rows % m_yNum;

	m_xPixelNum = (xRemainNum == 0 ? xPixelNum : xPixelNum + 1) << 1;
	m_yPixelNum = (yRemainNum == 0 ? yPixelNum : yPixelNum + 1) << 1;

	begin = (begin == -1 ? 0 : begin);
	end = (end == -1 ? m_yNum : end);

	for (int i = begin; i < end; i++)
	{
		CGrid *gridPtr = &m_imageGrid[i * m_xNum];

		for (int j = 0; j < m_xNum; j++)
		{
			gridPtr->x = xPixelNum * j + (j < xRemainNum ? j : xRemainNum);
			gridPtr->y = yPixelNum * i + (i < yRemainNum ? i : yRemainNum);

			gridPtr->width = (j < xRemainNum ? xPixelNum + 1 : xPixelNum);
			gridPtr->height = (i < yRemainNum ? yPixelNum + 1 : yPixelNum);

			gridPtr->size = gridPtr->width * gridPtr->height;
			gridPtr->label = i * m_xNum + j;


			float *meanVectorPtr = gridPtr->meanVector;
			uchar *pixelPtr = (uchar*)m_matData->data + (gridPtr->y * m_matData->cols + gridPtr->x) * m_matData->channels();

			for (int k = gridPtr->y; k < gridPtr->y + gridPtr->height; k++)
			{
				for (int l = gridPtr->x; l < gridPtr->x + gridPtr->width; l++)
				{
					meanVectorPtr[0] += ((float)l / (float)m_xPixelNum);
					meanVectorPtr[1] += ((float)k / (float)m_yPixelNum);
					meanVectorPtr[2] += (float)*pixelPtr++;
					meanVectorPtr[3] += (float)*pixelPtr++;
					meanVectorPtr[4] += (float)*pixelPtr++;
				}

				pixelPtr += (m_matData->cols - gridPtr->width) * m_matData->channels();
			}

			gridPtr->meanVector[0] /= (float)(gridPtr->size);
			gridPtr->meanVector[1] /= (float)(gridPtr->size);
			gridPtr->meanVector[2] /= ((float)(gridPtr->size) * NORMALIZE);
			gridPtr->meanVector[3] /= ((float)(gridPtr->size) * NORMALIZE);
			gridPtr->meanVector[4] /= ((float)(gridPtr->size) * NORMALIZE);

			gridPtr->updateMeanVector[0] = gridPtr->meanVector[0];
			gridPtr->updateMeanVector[1] = gridPtr->meanVector[1];
			gridPtr->updateMeanVector[2] = gridPtr->meanVector[2];
			gridPtr->updateMeanVector[3] = gridPtr->meanVector[3];
			gridPtr->updateMeanVector[4] = gridPtr->meanVector[4];

			gridPtr++;
		}
	}

	return true;
}

//  AssignLabel  //
void USEAQsuperpixel_TIP::AssignLabel(cv::Mat &labels, float omega, int begin, int end)
{
	static const int neighborNum = 9;

	static const int neighborIdxX[neighborNum] = { 0, 1, -1, 0, 0, -1, -1, 1, 1 };
	static const int neighborIdxY[neighborNum] = { 0, 0, 0, 1, -1, -1, 1, -1, 1 };

	CGrid *gridPtr;

	begin = (begin == -1 ? 0 : begin);
	end = (end == -1 ? m_yNum : end);

	for (int i = begin; i < end; i++)
	{
		gridPtr = &m_imageGrid[i * m_xNum];

		for (int j = 0; j < m_xNum; j++)
		{
			uchar *pixelPtr = (uchar*)m_matData->data + (gridPtr->y * m_matData->cols + gridPtr->x) * m_matData->channels();
			int *labelPtr = (int*)labels.data + (gridPtr->y * m_matData->cols + gridPtr->x);

			int gridSize = gridPtr->width * gridPtr->height;

			for (int k = gridPtr->y; k < gridPtr->y + gridPtr->height; k++)
			{
				for (int l = gridPtr->x; l <gridPtr->x + gridPtr->width; l++)
				{
					float maxGravitation = -FLT_MAX;
					int maxNeighborIdx = -1;

					float pixel[3];
					pixel[0] = (float)*pixelPtr++ / NORMALIZE;
					pixel[1] = (float)*pixelPtr++ / NORMALIZE;
					pixel[2] = (float)*pixelPtr++ / NORMALIZE;

					float xValue = (float)l / m_xPixelNum;
					float yValue = (float)k / m_yPixelNum;

					int meanVecIndex = 0;
					for (int m = 0; m < neighborNum; m++)
					{
						int neighborX = j + neighborIdxX[m];
						int neighborY = i + neighborIdxY[m];

						float colorTerm = FLT_MAX;
						float distanceTerm = FLT_MAX;
						if (neighborX >= 0 && neighborX < m_xNum && neighborY >= 0 && neighborY < m_yNum)
						{
							CGrid *neighborGridPtr = &m_imageGrid[neighborY * m_xNum + neighborX];
							float *neighborMeanPtr = neighborGridPtr->meanVector;

							int theNumOfMeanVec = 0;
							for (int dominColorInx = 0; dominColorInx < neighborGridPtr->dominateLabelCount; dominColorInx++)
							{
								float tempColorTerm = pow((pixel[0] - neighborGridPtr->dominateColorVector[dominColorInx][0]), 2) + pow((pixel[1] - neighborGridPtr->dominateColorVector[dominColorInx][1]), 2) + pow((pixel[2] - neighborGridPtr->dominateColorVector[dominColorInx][2]), 2);

								float tempDistanceTerm = pow(xValue - (neighborGridPtr->dominateSpatialVector[dominColorInx][0]), 2) + pow(yValue - (neighborGridPtr->dominateSpatialVector[dominColorInx][1]), 2);

								if (colorTerm  * (1.0f - (omega)) + (distanceTerm* omega)> tempColorTerm* (1.0f - (omega)) + (tempDistanceTerm * omega))
								{
									colorTerm = tempColorTerm;
									distanceTerm = tempDistanceTerm;
									theNumOfMeanVec = dominColorInx;
								}
							}

							float gravitation = 1 / (distanceTerm * omega + colorTerm * (1 - (omega)));

							if (maxGravitation < gravitation)
							{
								maxGravitation = gravitation;
								maxNeighborIdx = m;
								meanVecIndex = theNumOfMeanVec;
							}
						}
					}

					int maxNeighborX = j + neighborIdxX[maxNeighborIdx];
					int maxNeighborY = i + neighborIdxY[maxNeighborIdx];

					*labelPtr++ = (maxNeighborY * m_xNum + maxNeighborX) + (m_yNum * m_xNum)*meanVecIndex;

					if (maxNeighborIdx != 0)
					{
						CGrid *maxNeighborGridPtr = &m_imageGrid[maxNeighborY * m_xNum + maxNeighborX];
						float *maxMeanVectorPtr = maxNeighborGridPtr->updateMeanVector;
						float *meanVectorPtr = gridPtr->updateMeanVector;

						int newSize = maxNeighborGridPtr->size + 1;
						int size = gridPtr->size - 1;

						maxMeanVectorPtr[0] = (maxMeanVectorPtr[0] * maxNeighborGridPtr->size + xValue) / newSize;
						maxMeanVectorPtr[1] = (maxMeanVectorPtr[1] * maxNeighborGridPtr->size + yValue) / newSize;
						maxMeanVectorPtr[2] = (maxMeanVectorPtr[2] * maxNeighborGridPtr->size + pixel[0]) / newSize;
						maxMeanVectorPtr[3] = (maxMeanVectorPtr[3] * maxNeighborGridPtr->size + pixel[1]) / newSize;
						maxMeanVectorPtr[4] = (maxMeanVectorPtr[4] * maxNeighborGridPtr->size + pixel[2]) / newSize;

						meanVectorPtr[0] = (meanVectorPtr[0] * gridPtr->size - xValue) / size;
						meanVectorPtr[1] = (meanVectorPtr[1] * gridPtr->size - yValue) / size;
						meanVectorPtr[2] = (meanVectorPtr[2] * gridPtr->size - pixel[0]) / size;
						meanVectorPtr[3] = (meanVectorPtr[3] * gridPtr->size - pixel[1]) / size;
						meanVectorPtr[4] = (meanVectorPtr[4] * gridPtr->size - pixel[2]) / size;
						if (meanVectorPtr[0]<0)
							meanVectorPtr[0] = 0;
						if (meanVectorPtr[1] <0)
							meanVectorPtr[1] = 0;
						if (meanVectorPtr[2] <0)
							meanVectorPtr[2] = 0;
						if (meanVectorPtr[3] <0)
							meanVectorPtr[3] = 0;
						if (meanVectorPtr[4] <0)
							meanVectorPtr[4] = 0;

						maxNeighborGridPtr->updateSize = newSize;
						gridPtr->updateSize = size;
					}
				}
				pixelPtr += (m_matData->cols - gridPtr->width) * m_matData->channels();
				labelPtr += (labels.cols - gridPtr->width);
			}
			gridPtr++;
		}
	}
}


//  Cluster  //
int USEAQsuperpixel_TIP::Cluster(cv::Mat &matData, cv::Mat &labels, int spNum, float omega, bool tbbBoost)
{
	// initialize the number of superpixels
	double t;
	this->InitializeSp(matData, sqrt(spNum), sqrt(spNum));

	// perform color quantization & label assignment 
	if (tbbBoost) {
		cv::parallel_for_(cv::Range(0, m_yNum), BuildGridbyQuantizationParallel(this));
		cv::parallel_for_(cv::Range(0, m_yNum), AssignLabelParallel(m_labelData, omega, this));
	}
	else {
		this->BuildImageGridbyQuantization();
		this->AssignLabel(m_labelData, omega);
	}

	labels = m_labelData;

	int relableSpNum = m_xNum * m_yNum;
	cv::Mat m_labelsTemp;

	// perform label refinement
	this->LabelRefinement(m_labelData, m_labelsTemp, relableSpNum, m_xNum * m_yNum, omega);

	if (spData) delete[]spData;

	labels = m_labelsTemp;

	return relableSpNum;
}

int  USEAQsuperpixel_TIP::ColorQuntization(cv::Mat &matData, cv::Mat &labels, int qtzLV)
{
	int relableSpNum = 0;

	this->colorQuntization(matData, matData, labels, qtzLV);

	return relableSpNum;
}

void USEAQsuperpixel_TIP::colorQuntization(cv::Mat &_input, cv::Mat &_output, cv::Mat &_label, int dimension)
{
	if (_input.empty() == true)
		return;

	if ((_output.type() != CV_8UC3) || (_output.size() != _input.size()))
		_output = cv::Mat::zeros(_input.size(), _input.type());

	_label = cv::Mat::zeros(_input.size(), CV_32S);
	if (dimension == 0)
		return;

	int quntzNum = dimension;

	int qtzlv = dimension;
	int colorQtzize = ceil(255.0f / qtzlv);

	uchar *quntzColor = new uchar[quntzNum];
	for (int i = 0; i<quntzNum; i++)
	{
		quntzColor[i] = (256 / quntzNum) / 2 + ((256 / quntzNum)*i);
	}

	uchar *inputPtr = (uchar*)_input.data;
	uchar *outputtPtr = (uchar*)_output.data;
	int *_labelPtr = (int*)_label.data;
	for (int i = 0; i<_input.rows; i++)
	{
		for (int j = 0; j<_input.cols; j++)
		{
			int r, g, b;
			uchar c;

			r = (int)(*inputPtr++) / colorQtzize;
			c = quntzColor[r];
			(*outputtPtr++) = c;

			g = (int)(*inputPtr++) / colorQtzize;
			c = quntzColor[g];
			(*outputtPtr++) = c;

			b = (int)(*inputPtr++) / colorQtzize;
			c = quntzColor[b];
			(*outputtPtr++) = c;

			*_labelPtr++ = (int)(b*quntzNum*quntzNum) + (r*quntzNum) + g;
		}
	}
	if (quntzColor != nullptr)	delete[]quntzColor; quntzColor = nullptr;
}

//  LabelRefinement  //
bool USEAQsuperpixel_TIP::LabelRefinement(cv::Mat &m_labels, cv::Mat &m_nlabels, int &numlabels, const int &K, float omega)
{
	if (!m_nlabels.empty())
		m_nlabels.release();
	m_nlabels = cv::Mat(m_labelData.size(), m_labelData.type());
	int* labels = (int*)m_labels.data;
	int* nlabels = (int*)m_nlabels.data;

	const int dx4[4] = { -1, 0, 1, 0 };
	const int dy4[4] = { 0, -1, 0, 1 };

	int width = m_labels.cols;
	int height = m_labels.rows;
	const int sz = width * height;
	const int SUPSZ = sz / K;

	for (int i = 0; i < sz; i++) nlabels[i] = -1;
	int label(0);
	int* xvec = new int[sz];
	int* yvec = new int[sz];
	int oindex(0);
	int adjlabel(0);

	for (int j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			if (0 > nlabels[oindex])
			{
				nlabels[oindex] = label;
				xvec[0] = k;
				yvec[0] = j;

				int count(1);
				for (int c = 0; c < count; c++)
				{
					for (int n = 0; n < 4; n++)
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if ((x >= 0 && x < width) && (y >= 0 && y < height))
						{
							int nindex = y * width + x;

							if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}
					}
				}
				label++;
			}
			oindex++;
		}
	}


	numlabels = label;
	for (int i = 0; i < sz; i++)
	{
		labels[i] = nlabels[i];
		nlabels[i] = -1;
	}

	std::vector<int> drawCenterPoint;
	std::vector<int>traceRecordVector;

	spData = new SPdata[numlabels];
	int *labelIndex = new int[numlabels];
	for (int i = 0; i < numlabels; i++)  labelIndex[i] = i;

	label = 0;
	oindex = 0;

	for (int j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			if (0 > nlabels[oindex])
			{
				nlabels[oindex] = label;
				float pixelRGB[3] = { 0 };
				uchar *oindexPtr = m_matData->data + (oindex * 3);
				pixelRGB[0] = (float)*oindexPtr++;
				pixelRGB[1] = (float)*oindexPtr++;
				pixelRGB[2] = (float)*oindexPtr++;
				cv::Point loc(k, j);

				xvec[0] = k;
				yvec[0] = j;

				int count(1);
				for (int c = 0; c < count; c++)
				{
					for (int n = 0; n < 4; n++)
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if ((x >= 0 && x < width) && (y >= 0 && y < height))
						{
							int nindex = y * width + x;
							uchar *m_matDataPtr = m_matData->data + (nindex * 3);

							if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
							{
								pixelRGB[0] = pixelRGB[0] + (float)*m_matDataPtr++;
								pixelRGB[1] = pixelRGB[1] + (float)*m_matDataPtr++;
								pixelRGB[2] = pixelRGB[2] + (float)*m_matDataPtr++;
								loc.x = x + loc.x;
								loc.y = y + loc.y;
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}

							if (labels[oindex] != labels[nindex])
							{
								bool isFind = false;

								for (int i = 0; i<spData[label].neighborCount; i++)
								{
									if (spData[label].neighborLabel.get(i) == labels[nindex])
									{
										isFind = true;
										break;
									}

								}

								if (isFind != true)
								{
									spData[label].neighborLabel.add(labels[nindex]);
									spData[label].neighborCount++;

								}
							}
						}
					}
				}
				spData[label].label = &labelIndex[label];

				//normalized
				spData[label].meanColor[0] = (float)pixelRGB[0] / (float)count / NORMALIZE;
				spData[label].meanColor[1] = (float)pixelRGB[1] / (float)count / NORMALIZE;
				spData[label].meanColor[2] = (float)pixelRGB[2] / (float)count / NORMALIZE;
				spData[label].spCenter.x = (float)loc.x / (float)count / (float)width;
				spData[label].spCenter.y = (float)loc.y / (float)count / (float)height;

				spData[label].size = count;
				spData[label].belongThisLabel.add(label);

				spData[label].belongCount++;
				label++;
			}
			oindex++;
		}
	}


	bool  allDone = true;
	do {
		allDone = true;

		for (int p = 0; p < label; p++)
		{
			int i = *spData[p].label;
			if ((spData[i].size < SUPSZ * m_fRefinementMag) && spData[i].isKill == false)
			{
				allDone = false;
				float minColorDistance = FLT_MAX;
				int maxNeighborLabel = *spData[(int)spData[i].neighborLabel.get(0)].label;
				for (int t = 0; t < spData[i].neighborCount; t++)
				{
					int neighborLabel = *spData[(int)spData[i].neighborLabel.get(t)].label;
					if ((*spData[neighborLabel].label != *spData[i].label))
					{
						float colorDistance = pow((spData[i].meanColor[0] - spData[neighborLabel].meanColor[0]), 2)
							+ pow((spData[i].meanColor[1] - spData[neighborLabel].meanColor[1]), 2)
							+ pow((spData[i].meanColor[2] - spData[neighborLabel].meanColor[2]), 2);

						float spitalDistance = pow((spData[i].spCenter.x - spData[neighborLabel].spCenter.x), 2) + pow((spData[i].spCenter.y - spData[neighborLabel].spCenter.y), 2);

						float alpha = omega;
						if ((1.0 - alpha)*colorDistance + alpha * spitalDistance < minColorDistance)
						{
							minColorDistance = (1.0 - alpha)* colorDistance + alpha * spitalDistance;
							maxNeighborLabel = neighborLabel;
						}
					}
				}

				spData[maxNeighborLabel].meanColor[0] = ((spData[i].meanColor[0] * spData[i].size) + (spData[maxNeighborLabel].meanColor[0] * spData[maxNeighborLabel].size)) / (spData[maxNeighborLabel].size + spData[i].size);
				spData[maxNeighborLabel].meanColor[1] = ((spData[i].meanColor[1] * spData[i].size) + (spData[maxNeighborLabel].meanColor[1] * spData[maxNeighborLabel].size)) / (spData[maxNeighborLabel].size + spData[i].size);
				spData[maxNeighborLabel].meanColor[2] = ((spData[i].meanColor[2] * spData[i].size) + (spData[maxNeighborLabel].meanColor[2] * spData[maxNeighborLabel].size)) / (spData[maxNeighborLabel].size + spData[i].size);
				spData[maxNeighborLabel].size = spData[i].size + spData[maxNeighborLabel].size;

				for (int g = 0; g < spData[i].neighborCount; g++)
				{
					if (*spData[(int)spData[i].neighborLabel.get(g)].label != *spData[maxNeighborLabel].label)
					{
						spData[maxNeighborLabel].neighborLabel.add(spData[i].neighborLabel.get(g));
						spData[maxNeighborLabel].neighborCount++;
					}
				}
				for (int g = 0; g < spData[i].belongCount; g++)
				{
					*spData[(int)spData[i].belongThisLabel.get(g)].label = *spData[maxNeighborLabel].label;
					spData[maxNeighborLabel].belongThisLabel.add(spData[i].belongThisLabel.get(g));
					spData[maxNeighborLabel].belongCount++;
				}

				spData[maxNeighborLabel].belongCount;
				spData[i].isKill = true;
			}


		}
	} while (!allDone);

	for (int i = 0; i < sz; i++) nlabels[i] = -1;
	label = 0;
	oindex = 0;
	for (int j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			if (0 > nlabels[oindex])
			{
				nlabels[oindex] = *spData[label].label;

				xvec[0] = k;
				yvec[0] = j;

				int count(1);
				for (int c = 0; c < count; c++)
				{
					for (int n = 0; n < 4; n++)
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if ((x >= 0 && x < width) && (y >= 0 && y < height))
						{
							int nindex = y * width + x;

							if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = *spData[label].label;
								count++;
							}
						}

					}
				}

				label++;
			}
			oindex++;
		}
	}

	for (int i = 0; i < sz; i++)
	{
		labels[i] = nlabels[i];
		nlabels[i] = -1;
	}
	label = 0;
	oindex = 0;

	for (int j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			if (0 > nlabels[oindex])
			{
				nlabels[oindex] = label;

				xvec[0] = k;
				yvec[0] = j;
				int count(1);
				for (int c = 0; c < count; c++)
				{
					for (int n = 0; n < 4; n++)
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if ((x >= 0 && x < width) && (y >= 0 && y < height))
						{
							int nindex = y * width + x;

							if (0 > nlabels[nindex] && labels[oindex] == labels[nindex])
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;

							}
						}
					}
				}
				label++;
			}
			oindex++;
		}
	}

	for (int i = 0; i<sz; i++) labels[i] = nlabels[i];
	numlabels = label;

	if (xvec)
		delete[] xvec;
	if (yvec)
		delete[] yvec;
	if (labelIndex)
		delete[]labelIndex;
	return allDone;
}

void USEAQsuperpixel_TIP::LabelContourMask(cv::Mat &_labels, cv::Mat &result, int index, bool thick_line)
{
	cv::Mat _contour = cv::Mat(_labels.size(), CV_8UC3);
	if (_labels.type() != CV_32S)_labels.convertTo(_labels, CV_32S);
	if (result.empty() == true)result = cv::Mat::zeros(_labels.size(), CV_8UC3);
	if (result.type() != CV_8UC3) cv::cvtColor(result, result, CV_8UC3);

	int width = _labels.cols;
	int height = _labels.rows;

	const int dx8[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	int* labelPtr = (int*)_labels.data;
	uchar* contourPtr = (uchar*)_contour.data;
	uchar* resultPtr = (uchar*)result.data;
	for (int j = 0; j < height; j++)
	{
		for (int k = 0; k < width; k++)
		{
			int neighbors = 0;
			for (int i = 0; i < 8; i++)
			{
				int x = k + dx8[i];
				int y = j + dy8[i];

				if ((x >= 0 && x < width) && (y >= 0 && y < height))
				{
					int* neighborPtr = (int*)_labels.data + y * width + x;
					uchar* contourNeighborPtr = (uchar*)_contour.data + y * width + x;


					if (*labelPtr != *neighborPtr)
					{
						if (thick_line || !*contourNeighborPtr)
							neighbors++;
					}
				}
			}
			*labelPtr++;
			if (neighbors > 1)
			{
				if (index == 1)
				{
					(*contourPtr++) = -1;
					(*resultPtr++) = 0;
					(*resultPtr++) = 0;
					(*resultPtr++) = 0;
				}
				else
				{
					(*contourPtr++) = -1;
					(*resultPtr++) = 255;
					(*resultPtr++) = 255;
					(*resultPtr++) = 255;
				}

			}
			else
			{
				(*contourPtr++ = 0);
				(*resultPtr++);
				(*resultPtr++);
				(*resultPtr++);

			}


		}
	}
}

bool USEAQsuperpixel_TIP::Label2Color(cv::Mat &label, cv::Mat &output)
{
	if (label.empty()) return false;
	if (label.type() != CV_32S) { label.convertTo(label, CV_32S); }
	double min, max;
	cv::minMaxLoc(label, &min, &max);
	max++;
	if (max<0) { return false; }
	output = cv::Mat::zeros(label.size(), CV_8UC3);

	float start, stop;

	cv::vector<cv::Vec3b> color(max);
	cv::RNG rng = cv::theRNG();
	cv::Vec3b *colorPtr = color.data();
	for (int i = 0; i<max; i++)
	{
		cv::Vec3b newcolor(rng(255), rng(255), rng(255));
		(*colorPtr++) = newcolor;

	}

	cv::Vec3b* outputPtr = (cv::Vec3b*)output.data;
	int* labelPtr = (int*)label.data;
	for (int i = 0; i<output.rows; i++)
	{
		for (int j = 0; j<output.cols; j++)
		{
			if (*labelPtr >= 0 && *labelPtr<max)
			{
				cv::Vec3b c = color[*labelPtr];
				(*outputPtr++) = c;
			}
			else
			{
				(*outputPtr++) = 0;
			}
			*labelPtr++;
		}
	}
	color.clear();
	color.shrink_to_fit();
	return 0;
}

