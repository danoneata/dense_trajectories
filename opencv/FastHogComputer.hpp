inline
FastHogComputer::FastHogComputer(const IplImageWrapper& img, std::size_t nBins, bool fullOrientation, double pyramidScaleFactor, bool imgIs2DVector)
	: _integralHogBins(nBins), _fullOrientation(fullOrientation), _pyramidScaleFactor(pyramidScaleFactor)
{
	assert(img);
	assert(nBins > 1);
	if (!(img->nChannels == 1 || img->nChannels == 3) && !(img->nChannels == 2 && imgIs2DVector))
		throw std::runtime_error("FastHogComputer(): the image has to be either a gray, a color image, or a 2D vector field!");

	// construct the HOG integral images
	init(img, nBins, fullOrientation, pyramidScaleFactor, imgIs2DVector);
}

inline
FastHogComputer::FastHogComputer(const IplImagePyramid& imgPyr, std::size_t nBins, bool fullOrientation)
	: _integralHogBins(nBins), _fullOrientation(fullOrientation), _pyramidScaleFactor(imgPyr.getScaleFactor())
{
	assert(imgPyr);
	assert(nBins > 1);
	if (imgPyr.getImage((std::size_t)0)->nChannels != 1)
		throw std::runtime_error("FastHogComputer(): imgPyr is not a gray image!");
	if (imgPyr.getImage((std::size_t)0)->depth != IPL_DEPTH_64F)
		throw std::runtime_error("FastHogComputer(): format of imgPyr is not double!");

	// construct the HOG integral images
	init(imgPyr, nBins, fullOrientation);
}

inline
FastHogComputer::FastHogComputer(const FastHogComputer& fastHog)
	: _integralHogBins(fastHog._integralHogBins.begin(), fastHog._integralHogBins.end()),
	_fullOrientation(fastHog._fullOrientation), _pyramidScaleFactor(fastHog._pyramidScaleFactor)
{ }

inline
FastHogComputer::~FastHogComputer()
{ }

inline
FastHogComputer& FastHogComputer::operator=(const FastHogComputer& hogComputer) {
	_integralHogBins = hogComputer._integralHogBins;
	_fullOrientation = hogComputer._fullOrientation;
	_pyramidScaleFactor = hogComputer._pyramidScaleFactor;
	return *this;
}

inline
std::size_t FastHogComputer::numOfLevels() const {
	if (_integralHogBins.empty())
		return 0u;
	return _integralHogBins[0].numOfLevels();
}

inline
double FastHogComputer::getPyramidScaleFactor() const {
	return _pyramidScaleFactor;
}

inline
std::size_t FastHogComputer::numOfBins() const {
	return _integralHogBins.size();
}

inline
bool FastHogComputer::hasFullOrientation() const {
	return _fullOrientation;
}


