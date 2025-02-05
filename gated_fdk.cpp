/**
	CAVAREV 1.0 (CArdiac VAsculature Reconstruction EValuation) example reconstruction in C++.

	This is an example C++ application which is based on the data of the CAVAREV platform.
	It implements an ECG-gated FDK reconstruction as described on the website.
	It requires you to work with the post-processed projection data, the projection matrices
	and the heart phase data.
	
	The related data can be downloaded on our website: www.cavarev.com!

	Contact:
	Christopher Rohkohl
	Chair of Pattern Recognition, Department of Computer Science
	Friedrich-Alexander University Erlangen-Nuremberg
	Martensstr. 3, 91058 Erlangen, Germany.

	Email: contact@cavarev.com
**/



///-------------------------------------------------------------------------------------------
/// Basic includes and defines
///-------------------------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>
#include <cstring>
#include "cnpy.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif


template <class T> inline T round(float num)
{
	return static_cast<T>((num>0.0f) ? num+0.5f : num-0.5f);
}


///-------------------------------------------------------------------------------------------
/// Definition of basic structures and classes
///-------------------------------------------------------------------------------------------


/** Simple structure holding the information of a 3-D reconstruction */
struct Reco3D
{
	float *			volume;		///< memory for the 3-D reconstruction
	float			origin[3];	///< first voxel coordinate in world coordinate system [mm] -> isocenter: (0,0,0)
	unsigned int	size[3];	///< side lengths of the reconstruction in voxels
	float			voxelSize;	///< voxel side length in [mm]
};



/** A C++ class which handles loading of all relevant data and performs the ECG-gated FDK reconstruction.

	This class implements the core of the reconstruction routine. You should study this class extensively
	in order to understand how to deal with the CAVAREV data.
*/
class CavGatedFDK
{
	public:
		/** Create a reconstruction object with basic properties.
			\param N			number of projection images
			\param Sx			projection image width
			\param Sy			projection image height
			\param imageFile	path to the post-processed CAVAREV projection data
			\param matrixFile	path to the file containing the CAVAREV projection matrices
			\param phaseFile	path to the cardiac phase data. If no path is supplied, the heart phase is assumed 0 for all projection images.
		*/
		CavGatedFDK(unsigned int N, unsigned int Sx, unsigned int Sy, const std::string & imageFile, const std::string & matrixFile, const std::string & phaseFile="");
		~CavGatedFDK();


	public:
		/** Set the parameters of the cos^2 gating function.
			\param href		target reconstruction heart phase \in [0,1]
			\param width	The width of the weighting function \in (0,1].
		*/
		void setGating(float href=0.0f, float width=1.0f) { m_gating_href = href; m_gating_width=width; }


	public:
		/** Perform a 3-D reconstruction.
			\param reco		The parameters of the 3-D reconstruction.
		*/
		void doReconstruction(Reco3D & reco);


	protected:
		/** Perform 2-D interpolation with zero-border boundary condition in current projection image.
			\param x	zero-based column coordinate
			\param y	zero-based row coordinate
		*/
		float interpolate2D(float x, float y);

		/** Access a discrete image value with zero-border boundary condition in current projection image.
			\param i	zero-based column coordinate
			\param j	zero-based row coordinate
		*/
		inline float getImageValue(int i, int j);

		/** Evaluate the gating function for a certain heart phase.
			\param h	The heart phase for which to evaluate the gating function.
		*/
		float getGatingWeight(float h);


	private:
		unsigned int	m_N;				///< number of projection images
		unsigned int	m_Sx;				///< projection image width
		unsigned int	m_Sy;				///< projection image height
		float			m_gating_href;		///< reference heart phase for the gating function
		float			m_gating_width;		///< gating window width for the gating function
		float *			m_heartPhases;		///< array containing the heart phases for each projection image (all phases zero by default)
		float *			m_projMatrices;		///< array containing the projection matrices for each projection image
		float *			m_image;			///< current projection image
		std::ifstream	m_imageIfs;			///< file stream for reading the projection images
};



///-------------------------------------------------------------------------------------------
/// Implementation of the class CavGatedFDK
///-------------------------------------------------------------------------------------------

CavGatedFDK::CavGatedFDK(unsigned int N, unsigned int Sx, unsigned int Sy, const std::string & imageFile, const std::string & matrixFile, const std::string & phaseFile)
						: m_N(N), m_Sx(Sx), m_Sy(Sy)
{
	// allocate memory
	m_heartPhases	= new float[m_N];
	m_projMatrices	= new float[m_N*12];
	m_image			= new float[m_Sx*m_Sy];

	// initialize memory to zero
	memset(m_heartPhases, 0, m_N*sizeof(float));
	memset(m_projMatrices, 0, m_N*12*sizeof(float));
	memset(m_image, 0, m_Sx*m_Sy*sizeof(float));

	// read the phase data
	if (phaseFile.length()>0)
	{
		std::ifstream ifshp(phaseFile.c_str());

		for (unsigned int i=0; i<m_N; i++)
			ifshp >> m_heartPhases[i];

		ifshp.close();
	}

	// read the projection matrices
	std::ifstream ifspm(matrixFile.c_str(), std::ifstream::binary);
	ifspm.read((char*)(m_projMatrices), m_N*12*sizeof(float));
	ifspm.close();

	// open the input stream for the projection images
	m_imageIfs.open(imageFile.c_str(), std::ifstream::binary);

	// init with default gating parameters
	setGating();
}


CavGatedFDK::~CavGatedFDK()
{
	// cleanup memory
	delete [] m_heartPhases;
	delete [] m_projMatrices;
	delete [] m_image;

	// close streams
	m_imageIfs.close();
}


float CavGatedFDK::getGatingWeight(float h)
{
	float d = std::min(fabsf(h-m_gating_href), std::min(fabsf(h-m_gating_href+1.0f), fabsf(h-m_gating_href-1.0f)));

	if (d >= 0.5f*m_gating_width)
		return 0.0f;

	float cv = cosf(d / m_gating_width * (4.0f*atanf(1.0f)));
	return cv*cv;
}


void CavGatedFDK::doReconstruction(Reco3D & reco)
{
	// go to the beginning of the image file
	m_imageIfs.seekg(0, std::ifstream::beg);

	// initialize the reconstruction to zero
	unsigned int numelvol = reco.size[0]*reco.size[1]*reco.size[2];
	memset(reco.volume, 0, numelvol*sizeof(float));
	
	// for each projection image
	for (unsigned int i=0; i<m_N; i++)
	{
		std::cout << "Processing image " << (i+1) << " of " << m_N << std::endl;

		// load the post-processed projection image
		m_imageIfs.read((char*)(m_image), m_Sx*m_Sy*sizeof(float));
		// *****************************************************************
		// Save the projection image as a numpy array using cnpy
		std::vector<float> image_vector(m_Sx*m_Sy);
		for (unsigned int k=0; k<m_Sx*m_Sy; k++)
			{
				// std::cout << m_image[k] << std::endl;
				image_vector[k] = m_image[k];
				// if(m_image[k] != 0.0f)
				// {
				// 	std::cout << m_image[k] << std::endl;
				// }
			}
			// image_vector[k] = m_image[k];
			
		std::string filename = "projection_image_" + std::to_string(i) + ".npy";
		cnpy::npy_save(filename, image_vector.data(), {m_Sy,m_Sx}, "w");
		// *****************************************************************
		// determine the gating weight
		float lambda = getGatingWeight(m_heartPhases[i]);

		// skip the projection image if the weighting is zero
		if (lambda == 0.0f) continue;

		// multiply the projection image with the gating weight
		for (unsigned int k=0; k<m_Sx*m_Sy; k++)
			m_image[k] *= lambda;

		// determine the projection matrix
		float * A_i = m_projMatrices + i*12;

		// *****************************************************************
		// Save the projection image as a numpy array using cnpy
		std::vector<float> P(12);
		for (unsigned int k=0; k<12; k++)
			{
				// std::cout << m_image[k] << std::endl;
				P[k] = A_i[k];
				// if(m_image[k] != 0.0f)
				// {
				// 	std::cout << m_image[k] << std::endl;
				// }
			}
			// image_vector[k] = m_image[k];
			
		filename = "P" + std::to_string(i) + ".npy";
		cnpy::npy_save(filename, P.data(), {12}, "w");
		// *****************************************************************
		// perform the backprojection for each voxel x=(x0,x1,x2)
		float x[3];
		for (unsigned int ix2=0; ix2<reco.size[2]; ix2++)
		{
			x[2] = reco.origin[2] + (float)ix2*reco.voxelSize;
			for (unsigned int ix1=0; ix1<reco.size[1]; ix1++)
			{
				x[1] = reco.origin[1] + (float)ix1*reco.voxelSize;
				for (unsigned int ix0=0; ix0<reco.size[0]; ix0++)
				{
					x[0] = reco.origin[0] + (float)ix0*reco.voxelSize;

					// perform the backprojection of the voxel onto the
					// projection image
					float w_i =  A_i[8] * x[0] + A_i[9] * x[1] + A_i[10] * x[2] + A_i[11];
					float u_i = (A_i[0] * x[0] + A_i[1] * x[1] + A_i[2]  * x[2] + A_i[3]) / w_i;
					float v_i = (A_i[4] * x[0] + A_i[5] * x[1] + A_i[6]  * x[2] + A_i[7]) / w_i;

					// accumulate the voxel value
					reco.volume[ix2*reco.size[0]*reco.size[1] + ix1*reco.size[0] + ix0] += (float)(interpolate2D(u_i, v_i)/(w_i * w_i));
				}
			}
		}
	}
}

inline float CavGatedFDK::getImageValue(int i, int j)
{
	if (i>=0 && i<(int)m_Sx && j>=0 && j<(int)m_Sy)
		return m_image[j * m_Sx + i];

	return 0.0f;
}


float CavGatedFDK::interpolate2D(float x, float y)
{
	int i = (int)x;
	int j = (int)y;
	float alpha = x - (int)x;
	float beta  = y - (int)y;

	return (1.0f - alpha) * (1.0f - beta) * getImageValue(i  , j  )
	     +         alpha  * (1.0f - beta) * getImageValue(i+1, j  )
	     + (1.0f - alpha) *        beta   * getImageValue(i  , j+1)
	     +         alpha  *        beta   * getImageValue(i+1, j+1);
}



using namespace std;


///-------------------------------------------------------------------------------------------
/// Implementation of the main() program function.
/// It uses the previously defined reconstruction class.
/// It further calculates the reconstruction originl, sets the iso-center to the middle-most
/// volume voxel and saves the reconstruction in the CAVAREV 1.0 format.
///-------------------------------------------------------------------------------------------


// main function
int main(int argc, char **argv)
{
	// calling parameters
	// 0: reconstruction size in x dimension (number of matrix columns)
	// 1: reconstruction size in y dimension (number of matrix rows)
	// 2: reconstruction size in z dimension (number of slices)
	// 3: voxel size
	// 4: file path to post-processed projection data
	// 5: file path to projection matrices
	// 6: file path to heart phases
	// 7: gating window width in (0,1]
	// 8: reference heart phase in [0, 1]
	// 9: output file

	if (argc != 11)
	{
		cout << "Wrong calling syntax." << endl << endl << argv[0] << " [0] ... [9] " << endl
			 << "   0: reconstruction size in x dimension (number of matrix columns)" << endl
			 << "   1: reconstruction size in y dimension (number of matrix rows)" << endl
			 << "   2: reconstruction size in z dimension (number of slices)" << endl
			 << "   3: voxel size" << endl
			 << "   4: file path to post-processed projection data" << endl
			 << "   5: file path to projection matrices" << endl
			 << "   6: file path to heart phases" << endl
			 << "   7: gating window width in (0,1]" << endl
			 << "   8: reference heart phase in [0, 1]" << endl
			 << "   9: output file" << endl;

		return 1;
	}
	// a.out 512 512 512 0.5 ./cavarev_card.filtered.hq.bin ./cavarev.matrices.bin ./phasedata_card.txt 0.1 0.5 ./reco.vol
	// collect the input data from the command line


	// reconstruction properties
	Reco3D reco = {0};

	reco.size[0]	= atoi(argv[1]);
	reco.size[1]	= atoi(argv[2]);
	reco.size[2]	= atoi(argv[3]);
	reco.voxelSize	= (float)atof(argv[4]);
	reco.volume		= new float[reco.size[0]*reco.size[1]*reco.size[2]];

	// calculate the reconstruction coordinate of the first voxel (0,0,0) in
	// world coordinates. In this example we set the middle-most volume voxel
	// to the iso-center (origin of the world coordinate system (0,0,0)).
	reco.origin[0]	= -reco.voxelSize*0.5f*(float)(reco.size[0]-1); // x0 origin in world coordinates
	reco.origin[1]	= -reco.voxelSize*0.5f*(float)(reco.size[1]-1); // x1 origin in world coordinates
	reco.origin[2]	= -reco.voxelSize*0.5f*(float)(reco.size[2]-1); // x2 origin in world coordinates

	// file properties
	std::string file_projdata	= argv[5];
	std::string file_matrices	= argv[6];
	std::string file_phases		= argv[7];
	std::string file_output		= argv[10];

	// gating properties
	float gating_width = (float)atof(argv[8]);
	float gating_href  = (float)atof(argv[9]);

	// fixed data properties for cavarev
	const unsigned int N = 133; // number of projection images
	const unsigned int Sx = 960;
	const unsigned int Sy = 960;

	// perform the reconstruction
	CavGatedFDK gfdk(N, Sx, Sy, file_projdata, file_matrices, file_phases);

	gfdk.doReconstruction(reco);

	// save the reconstruction as 8-bit using min-max scaling
	float minval = std::numeric_limits<float>::max();
	float maxval = -std::numeric_limits<float>::max();
	unsigned int numelvol = reco.size[0]*reco.size[1]*reco.size[2];
	for (unsigned int i=0; i<numelvol; i++)
	{
		if (reco.volume[i] > maxval) maxval = reco.volume[i];
		if (reco.volume[i] < minval) minval = reco.volume[i];
	}

	unsigned char * volout = new unsigned char[numelvol];
	for (unsigned int i=0; i<numelvol; i++)
	{
		volout[i] = round<unsigned char>(255.0f * (reco.volume[i]-minval) / (maxval-minval));
	}

	// *****************************************************************
	// Save the reconstruction as a numpy array using cnpy
	std::vector<float> reco_vector(numelvol);
	for (unsigned int k=0; k<numelvol; k++)
	{
		reco_vector[k] = reco.volume[k];
	}
	std::string filename = "reconstruction.npy";
	cnpy::npy_save(filename, reco_vector.data(), {reco.size[2],reco.size[1],reco.size[0]}, "w");
	// *****************************************************************
	////// save the volume in the CAVAREV 1.0 format. All data is expected in little-endian format.
	// 3*float						3*4 bytes				first voxel x=(x0,x1,x2) origin in world coordinates
	// 3*unsigned int				3*4 bytes				reconstruction size (Sx0, Sx1, Sx2) in voxels
	// 1*float						1*4 bytes				voxel size in mm
	// Sx0*Sx1*Sx2*unsigned char	Sx0*Sx1*Sx2*1 bytes		reconstructed volume in row-major format
	std::ofstream ofsvol(file_output.c_str(), std::ofstream::binary);
	ofsvol.write((char*)(&(reco.origin[0])), 3*sizeof(float));
	ofsvol.write((char*)(&(reco.size[0])), 3*sizeof(unsigned int));
	ofsvol.write((char*)(&(reco.voxelSize)), sizeof(float));
	ofsvol.write((char*)(volout), numelvol*sizeof(unsigned char));
	ofsvol.close();

	// cleanup
	delete [] reco.volume;
	delete [] volout;

	return 0;
}

