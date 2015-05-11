#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "shadow.h"
#include "array2.h"

static float *global_d;

static int
compare_distances(const void *a,
                  const void *b)
{
    int i=*(const int*)a,
        j=*(const int*)b;
    float di=global_d[i],
          dj=global_d[j];
    if(di<dj) return -1;
    else if(di==dj) return 0;
    else return 1;
}

// t2 represents the normalized distance squared from the centre of the kernel; the value will be zero if t2>=1
static float
kernel(float t2)
{
    return sqr(1-t2);
}

static void
blot(Array2f &img,
     float x,
     float y,
     float r,
     float opacity)
{
    float r2=r*r;
    int ilo=(int)(x-r), ihi=(int)std::ceil(x+r),
        jlo=(int)(y-r), jhi=(int)std::ceil(y+r);
    if(ilo<0) ilo=0;
    if(ihi>img.ni-1) ihi=img.ni-1;
    if(jlo<0) jlo=0;
    if(jhi>img.nj-1) jhi=img.nj-1;
    for(int j=jlo; j<=jhi; ++j){
        float u=y-j;
        for(int i=ilo; i<=ihi; ++i){
            float v=x-i;
            float t=u*u+v*v;
            if(t<r2){
                img(i,j)*=std::exp(-opacity*kernel(t/r2));
            }
        }
    }
}

bool
compute_shadows(Array1<Vec3f> const &x,
                float radius,
                float opacity,
				Array1<float> const & dens,
                Vec3f const &light_position,
                float quality_factor,
                Array1<float> &illum,
                Array2f &shadow_map)
{
    const int n=x.n;
    illum.resize(n);
    if(n==0){
        return true; // nothing to do
    }

    // calculate distance from each point to light, also find bounding box for points
    Array1<float> d(n);
    Array1<int> order(n); // prepare order to render points in
    Vec3f xmin=x[0], xmax=x[0];
    for(int i=0; i<n; ++i){
        order[i]=i;
        update_minmax(x[i], xmin, xmax);
        d(i)=dist(x[i], light_position);
    }
    std::printf("  %d particles, bounding box (%f %f %f) - (%f %f %f)\n", n, xmin[0], xmin[1], xmin[2], xmax[0], xmax[1], xmax[2]);

    // get a rough bounding sphere (not optimal, but not awful)
    Vec3f centre=(xmin+xmax)/2;
    float r=dist(x[0], centre);
    for(int i=1; i<n; ++i){
        float ri=dist(x[i], centre);
        if(ri>r) r=ri;
    }
    std::printf("  bounding sphere centre at %f %f %f, radius %f\n", centre[0], centre[1], centre[2], r);

    // check light is adequately outside
    if(dist(light_position, centre)<1.25*r){
        std::printf("  light is too close to particles!!\n");
        return false;
    }

    // sort by distance from light
    global_d=d.data;
    qsort(order.data, n, sizeof(int), compare_distances);

    // figure out orthogonal rotation for viewing transform
    Vec3f camz=normalized(centre-light_position);
    Vec3f camy;
    if(std::fabs(camz[0])>=std::fabs(camz[1]) && std::fabs(camz[0])>=std::fabs(camz[2])){
        camy[0]=-camz[1];
        camy[1]=camz[0];
        camy[2]=0;
    }else{
        camy[0]=0;
        camy[1]=-camz[2];
        camy[2]=camz[1];
    }
    camy/=mag(camy);
    Vec3f camx=cross(camy,camz);

    // figure out desired image resolution
    int npix=int(quality_factor*r/radius);
    shadow_map.fill(npix, npix, 1.f);

    std::printf("  shadow map resolution %d x %d\n", npix, npix);
   
    // and scaling for perspective
    float scale=0.4*npix*dist(centre, light_position)/r;

    // render in order from closet to furthest, recording illumination as we go
    for(int i=0; i<n; ++i){
        int j=order[i];
        Vec3f p=x[j]-light_position;
        float u=dot(p, camx);
        float v=dot(p, camy);
        float w=dot(p, camz);
        assert(w>0);
        float px=scale*u/w+0.5f*npix,
              py=scale*v/w+0.5f*npix,
              pr=scale*radius/w;
        assert(px>=0 && px<npix);
        assert(py>=0 && py<npix);
        // look up illumination for this particle
        illum[j]=shadow_map((int)px, (int)py);
        // and then darken the image here
        blot(shadow_map, px, py, pr, opacity*dens[j]);
    }

    return true;
}

