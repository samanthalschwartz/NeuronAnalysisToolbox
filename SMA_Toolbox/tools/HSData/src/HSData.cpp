
#include "HSData.h"
#include <iomanip>
#include <sstream>
#include <cstdio>

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;

HSData::HSData(int _raw_sizeL, int _sizeY, int _sizeX, int _nT, float _gain,
               const ImageT &_bgMean, std::string _base_filename)
    : raw_sizeL(_raw_sizeL), sizeY(_sizeY), sizeX(_sizeX),
      nT(_nT), gain(_gain), bgMean(_bgMean),
      base_filename(_base_filename)
{
    sizeL=raw_sizeL-nBaselineRows;
    sensorPixels=sizeL*sizeY;
    raw_sensorPixels=raw_sizeL*sizeY;
    framePixels=sizeL*sizeY*sizeX;
    raw_framePixels=raw_sizeL*sizeY*sizeX;
    experimentPixels=sizeL*sizeY*sizeX*nT;
    raw_experimentPixels=raw_sizeL*sizeY*sizeX*nT;

    computeBaselineBgMean();
    omp_init_lock(&lock);
}
        

void HSData::loadAllFrames(PixelT *buf)
{
    #pragma omp parallel for
    for(int t=0;t<nT;t++) {
        HyperImageT im=makeHyperImage(buf+t*framePixels);
        loadFrame(t, im);
    }
}

void HSData::loadAllRawFrames(RawPixelT *buf)
{
    #pragma omp parallel for
    for(int t=0;t<nT;t++) {
        RawHyperImageT im= makeRawHyperImage(buf+t*raw_framePixels);
        loadRawFrame(t,im);
    }
}


void HSData::loadRawFrame(int t, RawHyperImageT &im)
{
    omp_set_lock(&lock);
    ICS *ics=openFile(t);
    IcsGetData(ics,im.memptr(), raw_framePixels*sizeof(RawPixelT));
    IcsClose(ics);
    omp_unset_lock(&lock);
}

void HSData::loadFrame(int t, HyperImageT &frame)
{
    RawPixelT rawData[raw_framePixels];
    omp_set_lock(&lock);
    ICS *ics=openFile(t);
    IcsGetData(ics, rawData, raw_framePixels*sizeof(RawPixelT));
    IcsClose(ics);
    omp_unset_lock(&lock);

    RawPixelT *rawSensorData=rawData; //rawSensorData points to the sizeL x sizeY column-major raw data
    PixelT *bgStart=bgMean.memptr(); //pointer to the bg data

    for(int x=0; x<sizeX;x++){
        PixelT baseline=computeBaseline(rawSensorData);
        PixelT *bg=bgStart; //rawSensorData points to the sizeL x sizeY column-major raw data
        for(int y=0; y<sizeY; y++) {
            for(int L=0; L<sizeL; L++) {
                //subtract bg; subtract baseline; correct gain.
                frame(L,y,x)=(rawSensorData[L]-bg[L]-baseline)*gain;
            }
            bg+=raw_sizeL; //Skip to next Y column
            rawSensorData+=raw_sizeL; //Skip to next Y column
        }
    }
}

std::string HSData::getFilename(int t) const
{
    ostringstream outs;
    outs << base_filename<<std::setw(5)<<std::setfill('0')<<t;
    return outs.str();
}

ICS* HSData::openFile(int t) const
{
    string filename=getFilename(t);
    ICS *ics;
    Ics_Error err=IcsOpen(&ics, filename.c_str(), "r");
    if (err!=IcsErr_Ok) {
        cout<<"ICS Error reading: \""<<filename<<"\"\n";
        cout<<"ICS Error Code:"<<err<<" Error Text:"<<IcsGetErrorText(err)<<endl;
        return nullptr;
    }
    return ics;
}


HSData::PixelT HSData::computeBaseline(RawPixelT *rawData) const
{
    assert(nBaselineRows==3);
    PixelT baseline=0.;
    int row1=raw_sizeL-1;
    int row2=raw_sizeL-2;
    int row3=raw_sizeL-3;
    for(int y=0; y<sizeY; y++) {
        baseline+=rawData[row1]+rawData[row2]+rawData[row3];//sum last 3 rows
        rawData+=raw_sizeL;//skip to next column
    }
    baseline/= nBaselineRows*sizeY; //Take the mean;
    baseline-=baselineBgMean; // Subtract the mean background in the last 3 rows
    return baseline;
}


void HSData::printFileInfo(ICS *ics, std::string &filename) const
{
    size_t data_size_bytes=IcsGetDataSize(ics);
    size_t data_size_pixels=IcsGetImageSize(ics);
    size_t pixel_size_bytes=IcsGetImelSize(ics);
    size_t dims[3];
    int nDims;
    Ics_DataType data_type;
    IcsGetLayout(ics, &data_type, &nDims, dims);
    cout<<"LibICS Version: "<<IcsGetLibVersion()<<"\n";
    cout<<"OMP thread: "<<omp_get_thread_num()<<"/"<<omp_get_num_procs()<<" (MAX: "<<omp_get_max_threads()<<")"<<"\n";
    cout<<"Opened file: "<<filename<<"\n";
    cout<<"  * Data size (bytes): "<<data_size_bytes<<"\n";
    cout<<"  * Data size (pixels): "<<data_size_pixels<<"\n";
    cout<<"  * Pixel size (bytes): "<<pixel_size_bytes<<"\n";
    cout<<"  * #Dims: "<<nDims<<"  [";
    for(int n=0;n<nDims;n++){
        cout<<dims[n];
        if(n<nDims-1) cout<<",";
    }
    cout<<"]\n";
    char coordinate_system[256];
    IcsGetCoordinateSystem(ics, coordinate_system);
    cout<<"  * Coordinate System: "<<coordinate_system<<"\n";
    for(int n=0;n<nDims;n++){
        char order[256], label[256];
        IcsGetOrder(ics, n, order, label);
        cout<<"  * Dim:"<<n<<" Order:"<<order<<" Label:"<<label<<"\n";
    }
}

void HSData::computeBaselineBgMean()
{
    //Want the mean of the last nBaselineRows of the sensor.
    baselineBgMean=arma::accu(bgMean.rows(raw_sizeL-nBaselineRows,raw_sizeL-1))/(nBaselineRows*sizeY);
}

std::ostream& operator<<(std::ostream &out, HSData &hsd)
{
    out<<"HSData:{rawSize(LYXT):["<<hsd.raw_sizeL<<","<<hsd.sizeY<<","<<hsd.sizeX<<","<<hsd.nT<<"}, "
       <<"size(LYXT):["<<hsd.sizeL<<","<<hsd.sizeY<<","<<hsd.sizeX<<","<<hsd.nT<<"}, "
       <<"#RawPixels{Sensor:"<<hsd.raw_sensorPixels<<", Frame:"<<hsd.raw_framePixels<<", Exp:"<<hsd.raw_experimentPixels<<"}, "
       <<"#Pixels{Sensor:"<<hsd.sensorPixels<<", Frame:"<<hsd.framePixels<<", Exp:"<<hsd.experimentPixels<<"}, "
       <<"gain:"<<hsd.gain<<"e-/ADU, mean(bgMean):"<<arma::accu(hsd.bgMean)/hsd.raw_sensorPixels<<"]"<<endl;
    return out;
}
