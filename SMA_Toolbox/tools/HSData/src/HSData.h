/** @file HSData.h
 * @author Mark J. Olah (mjo@@cs.unm.edu)
 * @date 07-17-2014
 * @brief This reads and gain corrects HS datasets in .ics format
 *
 *
 * .ics format.  From what I can tell the order of dimensions stored in the .ics file is [Y, L, X] but in row major format
 * Thus when we read in bytes in sequence to the column-major format of Armadillo/Matlab we end up with data indexed
 * as [L, Y, X].
 *
 */
#ifndef _HSDATA_H
#define _HSDATA_H

#include <cstdint>
#include <string>
#include <cassert>
#include <iostream>

#include <armadillo>
#include "libics.h"
#include "omp.h"

class HSData
{
public:
    typedef uint16_t RawPixelT;
    typedef float PixelT;
    typedef arma::Mat<PixelT> ImageT;
    typedef arma::Cube<PixelT> HyperImageT;
    typedef arma::Cube<RawPixelT> RawHyperImageT;
        
    const int nBaselineRows=3;
    int raw_sizeL; //This is the size before the extra baseline offset rows are removed
    int sizeL, sizeY, sizeX; //X is the scanned dimension
    int nT; //number of frames
    int sensorPixels; //sizeL*sizeY
    int raw_sensorPixels; //raw_sizeL*sizeY
    int framePixels;  //sizeL*sizeY*sizeX;
    int raw_framePixels;  //raw_sizeL*sizeY*sizeX;
    int experimentPixels; //sizeL*sizeY*sizeX*nT;
    int raw_experimentPixels; //raw_sizeL*sizeY*sizeX*nT;
    float gain; // photons/ADU
    ImageT bgMean; //sizeL x sizeY Normally: 128 x 64 Units ADU

    std::string base_filename; //Full directory plus filename except sequence number and .ics

    HSData(int _raw_sizeL, int _sizeY, int _sizeX, int _nT, float _gain,
            const ImageT &_bgMean, std::string _base_filename);

    /* Make new 3D data cubes for a frame data */
    RawHyperImageT makeRawHyperImage() const;
    HyperImageT makeHyperImage() const;
    RawHyperImageT makeRawHyperImage(RawPixelT *buf) const;
    HyperImageT makeHyperImage(PixelT *buf) const;

    /* Load all data into column-major 4D array buffers */
    void loadAllFrames(PixelT *buf);
    void loadAllRawFrames(RawPixelT *buf);

    /* Load data from frame t into column-major 3D Cube arrays */
    void loadFrame(int t, HyperImageT &frame);
    void loadRawFrame(int t, RawHyperImageT &im);

    friend std::ostream& operator<<(std::ostream &out, HSData &data);
private:
    float baselineBgMean; /* This is the mean of the BG image for the baseline correction rows */
    omp_lock_t lock; /* To serialize file ICS file reading which uses non-thread safe IO */

    /* Compute baseline corrections */
    void computeBaselineBgMean(); /* The BG Mean for the baseline rows is constant, so we precompute it */
    PixelT computeBaseline(RawPixelT *rawData) const;

    /* Manipulate filenames and open and read ICS files */
    std::string getFilename(int t) const;
    ICS* openFile(int t) const;
    void printFileInfo(ICS *ics, std::string &filename) const;
};


/**
 * @return An uninitialized raw_sizeL x sizeY x sizeX cube of RawPixelT
 * Raw images are indexed as [L, Y, X] in column major format.
 *
 */
inline
HSData::RawHyperImageT HSData::makeRawHyperImage() const
{
    return RawHyperImageT(raw_sizeL, sizeY, sizeX);
}

inline
HSData::HyperImageT HSData::makeHyperImage() const
{
    return HyperImageT(sizeL, sizeY, sizeX);
}

/**
 * @return An uninitialized raw_sizeL x sizeY x sizeX cube of RawPixelT
 * Raw images are indexed as [L, Y, X] in column major format.
 *
 */
inline
HSData::RawHyperImageT HSData::makeRawHyperImage(RawPixelT *buf) const
{
    assert(buf);
    return RawHyperImageT(buf, raw_sizeL, sizeY, sizeX, false, true);
}

inline
HSData::HyperImageT HSData::makeHyperImage(PixelT *buf) const
{
    assert(buf);
    return HyperImageT(buf, sizeL, sizeY, sizeX, false, true);
}

#endif /* _HSDATA_H */
