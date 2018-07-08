// USEQsuperpixel_TIP.cpp : 定義主控台應用程式的進入點。
//
#include "stdafx.h"
#include "USEAQsuperpixel_TIP.h"
#include <afxdialogex.h>
#include <fstream>
#include <iostream>
#include <direct.h>
#include <algorithm>
#include <functional>
#include<fstream>

#include "lazyGlobalVariable.h"

using namespace cv;
using namespace std;

bool FileFind(std::vector<std::string> &filePath, std::string setInputpath)
{
	CFileFind findFile;

	CString findPath = (setInputpath + "\\*.jpg").c_str();
	int nIsFind = findFile.FindFile((LPCTSTR)findPath);

	if (!nIsFind) {
		findPath = (setInputpath + "\\*.png").c_str();
		nIsFind = findFile.FindFile((LPCTSTR)findPath);
	}

	if (!nIsFind)
		return false;

	while (nIsFind)
	{
		nIsFind = findFile.FindNextFile();

		if (findFile.IsDots())
			continue;

		CStringA stra(findFile.GetFileName().GetBuffer(0)); //Convert CString to std::string
		filePath.push_back(stra.GetBuffer(0));
		stra.ReleaseBuffer();
	}

	return true;
}

struct Pred {
	Pred(int nth) : nth(nth) {};
	bool operator()(int k) { return k >= nth; }
	int nth;
};


int _tmain(int argc, _TCHAR* argv[])
{
	vector<string> findFileName;
	vector<float> timeRecoard;
	vector<int>	 vint_superpixelNum;

	//////////////////////////////////////////////////////////////////////////
	//set up parameters for superpixel extraction
	//m_fTheta - default as 4
	//m_fTheta is the number of each color channel be divided into.
	//m_fNumCandidates -  default as 10
	//m_fNumCandidates is maximum number of the adaptive sampling of superpixel candidates S in each spatially quantized region.
	//para_fOmega -  default as 0.01
	//para_fOmega controls the relative importance of the location information to the color information.
	//please refer to Section IV-B. Parameter Selection for more detail
	float m_fTheta = 4.0;
	float m_fNumCandidates = 10.0;
	float para_fOmega = 0.01;
	//////////////////////////////////////////////////////////////////////////	

	//////////////////////////////////////////////////////////////////////////
	//set up path, number of sp and output
	bool op_bStoreLabelImage = false;//whether to store the superpixel extraction result stored in int16
	bool op_bStoreColourImage = false;//whether to store the superpixel extraction result presented in different colour
	bool op_bStoreBoundedImage = false;//whether to store the superpixel extraction result presented by drawing the contour on the original image
	bool op_bShowResults = false;//whether to show results in different colour and contoured image on window

	//the path of image folder
	string path_strImageFolderPath = "E://TienJi//TIP_Dataset//BSDS//data//images//allimage";
	//string path_strImageFolderPath = "D://TienJi//2018_TIP//Dataset//Different Resolution//resolution_1440";
	//string path_strImageFolderPath = "E://TienJi//TIP_Dataset//SBD//images//test";

	//the path of result folder
	string path_strResultFolder = "E://TienJi//TIP_Dataset//Result//BSDS";

	// number of superpixels
	vint_superpixelNum.push_back(25);
	vint_superpixelNum.push_back(50);
	vint_superpixelNum.push_back(100);
	vint_superpixelNum.push_back(250);
	vint_superpixelNum.push_back(500);
	vint_superpixelNum.push_back(1000);
	vint_superpixelNum.push_back(2500);

	//scan jpgs or pngs in the dataset
	FileFind(findFileName, path_strImageFolderPath);

	//////////////////////////////////////////////////////////////////////////

	//int_imgs can be used to test if the exe works and the set up of path
	int int_imgs = findFileName.size();

	//whether to use parallel boosting
	bool bParallelBoosting = true;

	//create result folder
	_mkdir(path_strResultFolder.c_str());

	//set up bpf file which will be used in BSDS evaluation
	std::fstream txtFileTemp;
	txtFileTemp.open(path_strResultFolder + "//" + "bpf.txt", std::ios::out);
	//may require changes to run evaluation, the evaluation will scan the folder {BSDS500_root}\images\{mode}\*.jpg for imgs and {BSDS500_root}\groundTruth\{mode}\*.mat for gt
	txtFileTemp << "mode                allimage\n";
	txtFileTemp << "BSDS500_root " << path_strImageFolderPath << "\n";
	txtFileTemp << "nImages             500\n";
	txtFileTemp << "algResSavePath      " << path_strResultFolder << "\n";
	txtFileTemp << "algorithmCommand    [S time] = segment_box(I, pN)\n";
	for (int fl_intNumSPIdx = 0; fl_intNumSPIdx < vint_superpixelNum.size(); fl_intNumSPIdx++)
		txtFileTemp << "parameterset parametersetName USEAQ_" << vint_superpixelNum[fl_intNumSPIdx] << " pN " << vint_superpixelNum[fl_intNumSPIdx] << "\n";
	txtFileTemp.close();

	for (int fl_intNumSPIdx = 0; fl_intNumSPIdx < vint_superpixelNum.size(); fl_intNumSPIdx++)
	{
		// set the folder to store USEAQ superpixels
		string spNumberDir = path_strResultFolder + "//USEAQ_" + to_string(vint_superpixelNum[fl_intNumSPIdx]);
		_mkdir(spNumberDir.c_str());

		float time = 0.0f;

		// generate superpixels for all of the images in the folder
		for (int fl_intFileIdx = 0; fl_intFileIdx < int_imgs; fl_intFileIdx++)
		{
			USEAQsuperpixel_TIP doUSEAQ;

			string fileName = findFileName[fl_intFileIdx]; fileName.erase(fileName.size() - 4, 4);
			string resultPath = (spNumberDir + "//" + fileName);
			string labelPath, colorLabelPath, contourPath;
			Mat matInputImg = imread(path_strImageFolderPath + "//" + findFileName[fl_intFileIdx]);
			Mat matOriInput = matInputImg.clone();
			Mat matLabel, matColourImg, matBoundedImg;
			Mat matQtzLabel, matQtzColorLabel;

			//-------------------------------------------------------------------------------------------
			//In terms of colour, BuildImageGridbyQuantization only samples when over grid * (1.0f / float(m_fNumCandidates))
			doUSEAQ.m_fNumCandidates = (1.0f / float(m_fNumCandidates));
			//During the merge process in LabelRefinement, is has to be larger than grid * (1.0f / float(2.5)) to not merge with other superpixel
			doUSEAQ.m_fRefinementMag = (1.0f / float(2.5));
			doUSEAQ.m_fTheta = m_fTheta;
			cv::medianBlur(matInputImg, matInputImg, 3);

			double	t = (double)cv::getTickCount();
			doUSEAQ.Cluster(matInputImg, matLabel, vint_superpixelNum[fl_intNumSPIdx], para_fOmega, bParallelBoosting);

			//-------------------------------------------------------------------------------------------
			t = ((double)cv::getTickCount() - t) / cv::getTickFrequency();
			time = time + t;

			if (fl_intFileIdx % 500 == 0)
				std::cout << "Number of superpixel : " << vint_superpixelNum[fl_intNumSPIdx] << ", " << fl_intFileIdx << "___Totaltime: " << time << std::endl;

			if (op_bStoreLabelImage == true) {
				matLabel.convertTo(matLabel, CV_16U);

				labelPath = resultPath + ".png";
				imwrite(labelPath, matLabel);
			}

			if (op_bStoreBoundedImage) {
				matBoundedImg = matOriInput.clone();

				//draw the boundary of superpixels based on labels
				doUSEAQ.LabelContourMask(matLabel, matBoundedImg, -1);

				if (op_bShowResults == true)
					imshow("Contoured Image", matBoundedImg);

				contourPath = resultPath + "_bi.png";
				imwrite(contourPath, matBoundedImg);
			}

			if (op_bStoreColourImage) {
				doUSEAQ.Label2Color(matLabel, matColourImg);
				if (op_bShowResults == true)
					imshow("Colour Image", matColourImg);

				colorLabelPath = resultPath + "_cl.png";
				imwrite(colorLabelPath, matColourImg);
			}
			if (op_bShowResults == true)
				cvWaitKey(0);
		}
		cout << "------------Total time cost: " << time << "\t Average time: " << time / int_imgs << endl;
		timeRecoard.push_back(time / int_imgs);
	}
	std::fstream txtTimeRecord;
	txtTimeRecord.open(path_strResultFolder + "//" + "time.txt", std::ios::out);
	for (int w = 0; w < timeRecoard.size(); w++)
	{
		txtTimeRecord << "spNum:" << vint_superpixelNum[w] << "    aveTime:" << timeRecoard[w] << std::endl;
	}
	txtTimeRecord.close();
	timeRecoard.clear();

	system("pause");
	return 0;
}



