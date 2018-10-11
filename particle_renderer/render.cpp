#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "render.h"

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
blot(Array2<Vec3f> &image,
	 float x,
	 float y,
	 float r,
	 float opacity,
	 Vec3f const &rgb)
{
	float r2=r*r;
	int ilo=(int)(x-r), ihi=(int)std::ceil(x+r),
		jlo=(int)(y-r), jhi=(int)std::ceil(y+r);
	if(ilo<0) ilo=0;
	if(ihi>image.ni-1) ihi=image.ni-1;
	if(jlo<0) jlo=0;
	if(jhi>image.nj-1) jhi=image.nj-1;
	for(int j=jlo; j<=jhi; ++j){
		float u=y-j;
		for(int i=ilo; i<=ihi; ++i){
			float v=x-i;
			float t=u*u+v*v;
			if(t<r2){
				float a=std::exp(-opacity*kernel(t/r2));
				Vec3f color = a*Vec3f(image(i,j)[0],image(i,j)[1],image(i,j)[2])+(1-a)*rgb;
				image(i,j)=Vec3f(color[0],color[1],color[2]);
			}
		}
	}
}

// what is our shading model?
// - alpha compositing (atop) mediated by std::exp(-opacity*kernel(d/radius))
// - base colour + illum*light_rgb ? (we don't bother an angle-dependent scattering, or the correct 1/r^2 fall-off)

void
render_smoke(Array1<Vec3f> const &x,
			 Array1<float> const &illum,
			 Array1<Vec3f> const &base_rgb, // particle colour with just ambient illumination
			 Array1<Vec3f> const &lit_rgb, // particle colour in full lighting
			 float radius,
			 float opacity,
			 Array1<float> const &den, 
			 Vec3f const &light_position,
			 Vec3f const &background_rgb,
			 Vec3f const &camera_position,
			 Vec3f const &camera_target,
			 float const camera_perspective,
			 Array2<Vec3f> &image)
{
	//Array2<Vec4f> image2;
	//image2.resize(image.ni,image.nj,Vec4f(0,0,0,0));
	const int n=x.n;
	image.assign(image.ni, image.nj, background_rgb);

	// sort particles by distance from camera
	Array1<float> d(n);
	Array1<int> order(n); // prepare order to render points in
	for(int i=0; i<n; ++i){
		order[i]=i;
		d(i)=dist(x[i], camera_position);
	}
	global_d=d.data;
	qsort(order.data, n, sizeof(int), compare_distances);

	// figure out camera transform stuff
	Vec3f camz=normalized(camera_position-camera_target);
	Vec3f camx=normalized(cross(Vec3f(0,0,1),camz)); 
	Vec3f camy=cross(camz, camx);

	// render back to front
	for(int i=n-1; i>=0; --i){
		int j=order[i];
		Vec3f p=x[j]-camera_position;
		float u=dot(p, camx),
			v=dot(p, camy),
			w=dot(p, camz);
		if(w>=0) continue; // behind the camera
		float px=(-camera_perspective*u/w)*image.nj+0.5f*image.ni,
			py=(-camera_perspective*v/w+0.5f)*image.nj,
			pr=(-camera_perspective*radius/w)*image.nj;
		blot(image, px, py, pr, opacity*den[j], (1-illum[j])*base_rgb[j]+illum[j]*lit_rgb[j]);
	}
	/*float max_illum = 0;
	for (int i=0; i<image.size();i++)
	{
		Vec3f color = Vec3f(image2[i][0],image2[i][1],image2[i][2]);
		max_illum = max(mag(color),max_illum);
	}*/
	//for (int i=0; i<image.size();i++)
	//{
	//	//float alpha = image2[i][3]/max_alpha;
	//	//alpha = min(alpha,1.0f);
	//	
	//	Vec3f color = Vec3f(image2[i][0],image2[i][1],image2[i][2]);
	//	float alpha = mag(color)/max_illum;
	//	image[i] = color;
	//}
	//
}
