#pragma once
#include <afxtempl.h>
#include <opencv\cv.h>
#include <opencv\highgui.h>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>

#include "ArrayList.h"
#include <fstream>
#include <iostream>

class USEAQsuperpixel_TIP
{
private:
	struct CGrid
	{
		int x;
		int y;
		int width;
		int height;
		int label;
		float meanVector[8];
		float dominateColorVector[10][3];
		float dominateSpatialVector[10][2];
		float dominateLabelSize[10];
		float dominateLabelIdx[10][3];
		float updateMeanVector[8];
		int dominateLabelCount;
		int size;
		int updateSize;
		CGrid() : meanVector(), updateMeanVector(), dominateLabelCount(1) {}
		~CGrid() {}
	};


public:
	class SPdata
	{
	public:
		cv::Point2f spCenter;
		int neighborCount;
		ArrayList neighborLabel;
		ArrayList belongThisLabel;
		int belongCount;
		float meanColor[3];
		int *label;
		int size;
		bool isKill;
		SPdata() :neighborLabel(500), belongThisLabel(100), neighborCount(0), belongCount(0), isKill(false) {}
		~SPdata() {}
	};

public:
	int dominateColorNum;
	float m_fRefinementMag;
	float m_fNumCandidates;
	int  m_fTheta;

private:
	cv::Mat *m_matData;
	cv::Mat m_inputTemp;
	cv::Mat m_labelData;
	cv::Mat m_quanzitationMap;
	int m_xNum;
	int m_yNum;
	int m_xPixelNum;
	int m_yPixelNum;
	CGrid *m_imageGrid;
	SPdata *spData;

public:
	USEAQsuperpixel_TIP(void);
	~USEAQsuperpixel_TIP(void);

private:
	class BuildGridParallel : public cv::ParallelLoopBody
	{
	private:
		USEAQsuperpixel_TIP * _KDClusterNspPtr;

	public:
		BuildGridParallel(USEAQsuperpixel_TIP *KDClusterNspPtr)
			:_KDClusterNspPtr(KDClusterNspPtr) {}

		void operator() (const cv::Range& range) const
		{
			_KDClusterNspPtr->BuildImageGrid(range.start, range.end);
		}
	};

	class BuildGridbyQuantizationParallel : public cv::ParallelLoopBody
	{
	private:
		USEAQsuperpixel_TIP * _KDClusterNspPtr;

	public:
		BuildGridbyQuantizationParallel(USEAQsuperpixel_TIP *KDClusterNspPtr)
			:_KDClusterNspPtr(KDClusterNspPtr) {}

		void operator() (const cv::Range& range) const
		{
			_KDClusterNspPtr->BuildImageGridbyQuantization(range.start, range.end);
		}
	};

	class AssignLabelParallel : public cv::ParallelLoopBody
	{
	private:
		cv::Mat _labels;
		float _omega;

		USEAQsuperpixel_TIP *_KDClusterNspPtr;

	public:
		AssignLabelParallel(cv::Mat& labels, double omega, USEAQsuperpixel_TIP *KDClusterNspPtr)
			:_labels(labels), _omega(omega), _KDClusterNspPtr(KDClusterNspPtr) {}

		void operator() (const cv::Range& range) const
		{
			cv::Mat dstStripe = _labels;
			_KDClusterNspPtr->AssignLabel(dstStripe, _omega, range.start, range.end);
		}
	};

private:
	bool InitializeSp(cv::Mat &matData, int xNum, int yNum);
	bool BuildImageGrid(int begin = -1, int end = -1);
	bool BuildImageGridbyQuantization(int begin = -1, int end = -1);
	void colorQuntization(cv::Mat &_input, cv::Mat &_output, cv::Mat &_label, int dimension);
	void AssignLabel(cv::Mat &labels, float omega, int begin = -1, int end = -1);
	bool LabelRefinement(cv::Mat &m_labels, cv::Mat &m_nlabels, int &numlabels, const int &K, float omega);

public:
	int Cluster(cv::Mat &matData, cv::Mat &labels, int spNum, float omega = 0.01, bool tbbBoost = true);
	int  ColorQuntization(cv::Mat &matData, cv::Mat &labels, int qtzLV = 3);
	bool Label2Color(cv::Mat &label, cv::Mat &output);
	void LabelContourMask(cv::Mat &_labels, cv::Mat &result, int index, bool thick_line = true);
};

