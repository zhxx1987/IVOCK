#include "bfstream.h"
#include "write_sgi.h"

// SGI file IO ===============================================================

bool write_sgi(const Array2f &img, bool high_precision, const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);
    bofstream output;
    output.vopen(filename_format, ap);
    va_end(ap);
    if(!output.good()) return false;

    assert(sizeof(short)==2 && sizeof(int)==4);
    output.set_big_endian();
    output<<(short)474<<(char)0<<(char)(high_precision?2:1);
    output<<(unsigned short)2<<(unsigned short)img.ni<<(unsigned short)img.nj<<(unsigned short)1;
    output<<0<<(high_precision?65535:255);
    output.write_zero(492);

    // now write the actual data
    for(int j=0; j<img.nj; ++j){
        for(int i=0; i<img.ni; ++i){
            float scaled_value=img(i,j)*(high_precision?65536:256);
            if(scaled_value<0) scaled_value=0;
            else if(!high_precision && scaled_value>255) scaled_value=255;
            else if(high_precision && scaled_value>65535) scaled_value=65535;
            if(!output.good()) return false;
            if(high_precision) output<<(unsigned short)scaled_value;
            else output<<(unsigned char)scaled_value;
        }
    }
    return !output.fail();
}

bool write_sgi(const Array2<Vec3f> &img, bool high_precision, const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);
    bofstream output;
    output.vopen(filename_format, ap);
    va_end(ap);
    if(!output.good()) return false;

    assert(sizeof(short)==2 && sizeof(int)==4);
    output.set_big_endian();
    output<<(short)474<<(char)0<<(char)(high_precision?2:1);
    output<<(unsigned short)3<<(unsigned short)img.ni<<(unsigned short)img.nj<<(unsigned short)3;
    output<<0<<(high_precision?65535:255);
    output.write_zero(492);

    // now write the actual data
    for(int channel=0; channel<3; ++channel){
        for(int j=0; j<img.nj; ++j){
            for(int i=0; i<img.ni; ++i){
                float scaled_value=img(i,j)[channel]*(high_precision?65536:256);
                if(scaled_value<0) scaled_value=0;
                else if(!high_precision && scaled_value>255) scaled_value=255;
                else if(high_precision && scaled_value>65535) scaled_value=65535;
                if(!output.good()) return false;
                if(high_precision) output<<(unsigned short)scaled_value;
                else output<<(unsigned char)scaled_value;
            }
        }
    }
    return !output.fail();
}
