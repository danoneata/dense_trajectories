#ifndef IPLIMAGEPYRAMID_H_
#define IPLIMAGEPYRAMID_H_

// STL
#include <vector>

// my stuff
#include "IplImageWrapper.h"


class IplImagePyramid {
protected:
	std::vector<IplImageWrapper> _imagePyramid;
	std::vector<double> _scaleFactors;
	std::vector<double> _scaleFactorsInv;
	std::vector<double> _xScaleFactors;
	std::vector<double> _xScaleFactorsInv;
	std::vector<double> _yScaleFactors;
	std::vector<double> _yScaleFactorsInv;
	double _scaleFactor;
	double _epsilon;

protected:
//	IplImagePyramid(std::vector<IplImageWrapper> imagePyramid, std::vector<double> correctScaleFactors);

public:
	IplImagePyramid();

	IplImagePyramid(const IplImagePyramid& pyramid);

	/**
	 * NOTE: the image is referenced on the lowest pyramid level!
	 */
	IplImagePyramid(IplImageWrapper image, double scaleFactor);

	/**
	 * Build an empty pyramid (pixel values are set to zero).
	 */
	IplImagePyramid(CvSize initSize, int depth, int nChannels, double scaleFactor);

	~IplImagePyramid();

	IplImagePyramid& operator=(const IplImagePyramid& pyramid);

	operator const bool() const;

	operator bool();

	std::size_t numOfLevels() const;

	double getScaleFactor() const;

	/**
	 * round == 0  =>  take the closest level (i.e., rounding)
	 * round < 0   =>  take the next level with a smaller factor (i.e., flooring)
	 * round > 0   =>  take the next level with a bigger factor (i.e., ceiling)
	 */
	std::size_t getIndex(double scaleFactor, int round = 0) const;

	double getScaleFactor(std::size_t index) const;

	double getScaleFactorInv(std::size_t index) const;

	double getXScaleFactor(std::size_t index) const;

	double getXScaleFactorInv(std::size_t index) const;

	double getYScaleFactor(std::size_t index) const;

	double getYScaleFactorInv(std::size_t index) const;

	IplImageWrapper& getImage(std::size_t index);

	const IplImageWrapper& getImage(std::size_t index) const;

	/**
	 * @param round  see getIndex()
	 */
	IplImageWrapper getImage(double scaleFactor, int round = 0);
	const IplImageWrapper& getImage(double scaleFactor, int round = 0) const;

	/**
	 * rebuilds the pyramid (re-using the already allocated space) with the given 
	 * image
	 * NOTE: this image needs to have the exact sames size as the initial scale
	 */
	void rebuild(IplImageWrapper image);

private:
	void init(IplImageWrapper image, double scaleFactor);

	/**
	 * Build an empty pyramid (pixel values are set to zero).
	 */
	void init(CvSize initSize, int depth, int nChannels, double scaleFactor);

};

#include "IplImagePyramid.hpp"

#endif /*IPLIMAGEPYRAMID_H_*/
