#include "HSData.h"

using std::cout;
using std::endl;

int raw_sizeL=128;
int sizeY=64;
int sizeX=128;
int nT=250;
float gain=0.064834386110306; // e-/ADU
HSData::ImageT bgMean(raw_sizeL, sizeY, arma::fill::ones);

//std::string base_filename("/massive/share/LidkeLab/Data/M_2012_12_13A/HEK_EP4-HA_Cell04_100nM_PGE2_06min-2012-12-13-16-28-39_");
// std::string base_filename("/home/mjo/LidkeLab/Data/M_2012_12_13/M_2012_12_13A/HEK_EP4-HA_Cell04_100nM_PGE2_06min-2012-12-13-16-28-39_");
std::string base_filename("/home/mjo/LidkeLab/Data/OceanNanoBeads/OceanNano_530_575_600_665_try2-2013-12-03-17-59-52_");

int main()
{
    bgMean*=96.;
    HSData hsd(raw_sizeL, sizeY, sizeX, nT, gain, bgMean, base_filename);

    HSData::RawPixelT *raw_pixels=new HSData::RawPixelT[hsd.raw_experimentPixels];
    hsd.loadAllRawFrames(raw_pixels);
    float raw_sum=0.0;
    for(long n=0; n<hsd.raw_experimentPixels; n++) {
       raw_sum+=raw_pixels[n];
    }
    cout<<"Got sum: "<<raw_sum<<"\n";
    delete[] raw_pixels;

    cout<<"HSD: "<<hsd;
    HSData::PixelT *pixels=new HSData::PixelT[hsd.experimentPixels];
    hsd.loadAllFrames(pixels);


    float sum=0.0;
    for(long n=0; n<hsd.experimentPixels; n++) {
        sum+=pixels[n];
    }
    cout<<"Got sum: "<<sum<<"\n";
    delete[] pixels;

    return 1;
}







