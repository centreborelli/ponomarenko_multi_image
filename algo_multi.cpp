/*
 *
 * This file is part of the Ponomarenko Noise Estimation algorithm.
 *
 * Copyright(c) 2011 Miguel Colom.
 * miguel.colom@cmla.ens-cachan.fr
 *
 * This file may be licensed under the terms of of the
 * GNU General Public License Version 2 (the ``GPL'').
 *
 * Software distributed under the License is distributed
 * on an ``AS IS'' basis, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied. See the GPL for the specific language
 * governing rights and limitations.
 *
 * You should have received a copy of the GPL along with this
 * program. If not, go to http://www.gnu.org/licenses/gpl.html
 * or write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include <iostream>
#include <assert.h>
#include <string>
#include <dirent.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "algo.h"
#include "curve_filter.h"
//
#include "framework/CFramework.h"
#include "framework/CImage.h"
#include "framework/libparser.h"
#include "framework/operations.cpp"
#include "framework/CHistogram.cpp"


using namespace std;

//! A data structure to store the global id of a block within a sequence
// of images
struct BlockId
{
  int image;
  int block_origin;
  int bin;
};

template <typename T>
void print_arr(const T* data, int len)
{
  for (size_t i = 0; i < len; ++i) {
    cout << data[i] << ", ";
  }
  cout << endl;
}


//! Computes the delta matrix and returns the normalization factor theta
/*!
  \param *delta Delta matrix (mask for the low/high freqs in the block)
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \return theta Normalization factor for the matrix delta
*/
int compute_delta(float *delta, int w, int T)
{
  int theta = 0;
  for (int j = 0; j < w; j++)
    for (int i = 0; i < w; i++)
    {
      int value = (i + j < T && i + j != 0 ? 1 : 0);
      delta[j * w + i] = value;
      theta += value;
    }
  return theta;
}

//! Computes the set of variances computed form the low-frequency coefficients of the given blocks
/*!
  \param *VL Output set of variances
  \param M number of blocks taken into account
  \param w Block side
  \param *delta Delta matrix (mask for the low/high freqs in the block)s
  \param **blocks_ptr List of pointers to the blocks
  \param theta Normalization factor for the matrix delta
*/
void compute_VL(float *VL, int M, int w, float *delta, float **blocks_ptr,
                int theta)
{
  for (int m = 0; m < M; m++)
  {
    float *block = blocks_ptr[m];
    VL[m] = 0;

    for (int j = 0; j < w; j++)
    {
      for (int i = 0; i < w; i++)
        if (delta[j * w + i] != 0)
          VL[m] += pow(block[j * w + i], 2);
    }
    VL[m] /= theta;
  }
}

//! Computes the set of variances computed from the high-frequency coefficients of the given blocks
/*!
  \param *VH Output set of variances
  \param **blocks_ptr List of pointers to the blocks
  \param *indices_VL Sorting indices for the blocks_ptr list (by low-freqs)
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \param K Number of blocks that should be used
  \return Length of the returned variances list
*/
int compute_VH(float *VH, float **blocks_ptr, int *indices_VL, int w,
               int T, int K)
{
  int VH_count = 0;

  //#pragma omp parallel for
  for (int q = 0; q < w * w; q++)
  {
    int j = q / w;
    int i = q - j * w;

    if (i + j >= T)
    {
      float s = 0.0;
      for (int k = 0; k < K; k++)
      {
        float *block = blocks_ptr[indices_VL[k]];
        s += pow(block[q], 2); // q == j*w+i
      }
      VH[VH_count++] = s / K;
    }
  }
  return VH_count;
}

//! Computes the optimal K parameter using Ponomarenko's original article loop
/*!
  \param M Number of variance values in VL to use
  \param *VL List of variances obtained for low-freq coefficients
  \return The optimal K
*/
int get_optimal_K_ponom_orig(int M, float *VL)
{
  int K = sqrt(M);
  //
  for (int i = 0; i < 7; i++)
  {
    float U = 1.3 * VL[K / 2];

    int m_min = arg_find<float>(U, VL, M);

    int K1 = m_min;
    if (K1 > 0)
      K = K1;
  }

  // Set K = K / 5 to provide robustness
  int K1 = int(K / 5.0);
  if (K1 > 0)
    K = K1;

  return K;
}

//! Return the optimal T parameter according to the given block side
/*!
  \param w Block side
  \return The optimal T parameter
*/
int get_T(int w)
{
  switch (w)
  {
  case 3:
    return 3;
  case 4:
    return 3;
  case 5:
    return 5;
  case 7:
    return 8;
  case 8:
    return 9;
  case 11:
    return 13;
  case 15:
    return 17;
  case 21:
    return 24;
  default:
    PRINT_ERROR("Unknown block side: %d\n", w);
    exit(-1);
  }
}



//! Build an index map from valid blocks to raw blocks of an image,
//! so that the origin of the i-th block in a valid block list corresponds to the valid_coords[i]-th 
//! pixel in the raw image
/*!
  \param *mask Input valid mask of an image, 0 for valid, 1 for invalid
  \param **valid_coords pointer to output index map
  \param Nx Length of a row in the image
  \param Ny Length of a column in the image
  \param w Block side
  \param num_blocks Number of valid blocks
*/
void make_valid_coords(const int* mask, unsigned* valid_coords,
                       int Nx, int Ny) {
  int count_coords = 0;
  //    
  for (int i = 0; i < Nx*Ny; i++) {
    if (mask[i] == 0) {
      valid_coords[count_coords++] = i;

    }
  }
}


//! Reads all valid blocks (all neighbor pixels are different when the mask
//! is active) in the image
/*!
  \param *D Output list of blocks
  \param *u Input image
  \param Nx Length of a row in the image
  \param Ny Length of a column in the image
  \param w Block side
  \param num_blocks Number of blocks
  \param *mask Mask to determine if a pixel is valid or not
  \return Number of valid block copied into the output list
*/
void read_all_valid_blocks(float *D,
                           float *u,
                           int Nx, int Ny,
                           int w, unsigned num_blocks, int *mask) {  
  if (mask == NULL) {
    const int w2 = w * w;
    int q = 0;
    for (int y = 0; y < Ny - w + 1; ++y) {
      for (int x = 0; x < Nx - w + 1; ++x) {
        for (int j = 0; j < w; ++j) {
          for (int i = 0; i < w; ++i) {
            D[q*w2+j*w+i] = u[(j+y)*Nx+i+x];
          }
        }
        ++q;
      }
    }
  }
  
  else {
    unsigned* valid_coords = new unsigned[num_blocks];
    make_valid_coords(mask, valid_coords, Nx, Ny);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (unsigned q = 0; q < num_blocks; q++) {
      int addr = valid_coords[q];      

      for (int j = 0; j < w; j++) {
        for (int i = 0; i < w; i++) {
            D[q*w*w+j*w+i] = u[j*Nx+i+addr];
        }
      }
    }
    delete[] valid_coords;
  }
}


//! Computes the mean of all given blocks
/*!
  \param *means Output list of means of blocks
  \param *blocks Input list of blocks to compute their means
  \param w Block side
  \param num_blocks Number of block in the input list
*/
void compute_means(float *means, float *blocks, int w, int num_blocks)
{
  const float ONE_DIV_w2 = 1.0 / (w * w);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int b = 0; b < num_blocks; b++)
  {
    float mean = 0.0;
    for (int p = 0; p < w * w; p++)
    {
      mean += blocks[b * w * w + p];
    }
    mean *= ONE_DIV_w2;
    means[b] = mean;
  }
}

int get_max(int *data, int N)
{
  int max = data[0];
  for (int i = 1; i < N; i++)
    if (data[i] > max)
      max = data[i];
  return max;
}

void copy(float *dest, float *orig, int N)
{
  for (int i = 0; i < N; i++)
    dest[i] = orig[i];
}

//! Determines if the given string corresponds to the custom percentile code
/*!
  \param *num_str Input string
  \return true if the input string corresponds to the custom percentile code or false if not.
*/
bool is_custom_percentile(const char *num_str)
{
  char buffer[1024];
  float value = atof(num_str);
  sprintf(buffer, "%.4f", value);
  return strcmp(buffer, "0.0000") != 0;
}

//! Returns the mean of the data associated to the top K indices
/*!
  \param mean_method Method to compute the mean (1: mean of means, 2: median of means)
  \param K Number of selected elements for computing the mean
  \param indices List of indices associated to the data
  \param data raw data
  \return The mean of the bin
*/
float get_top_K_mean(int mean_method, int K, int *indices,
                     float *data)
{
  float mean;

  float *values = new float[K + 1];
  for (int i = 0; i <= K; i++)
    values[i] = data[indices[i]];

  switch (mean_method)
  {
  case 1:
  { // mean of means
    mean = 0.0;
    for (int i = 0; i <= K; i++)
    {
      mean += values[i];
    }
    mean /= (K + 1);
    break;
  }
  case 2:
  { // median of means
    mean = median<float>(values, K + 1);
    break;
  }
  default:
  {
    PRINT_ERROR("Unknown mean method: %d\n", mean_method);
    exit(-1);
  }
  }
  delete[] values;

  return mean;
}

//! In-place Normalization of the FFTW output in order to get a orthonormal 2D DCT-II
/*!
  \param *blocks Input/output list of transformed blocks
  \param w Block side
  \param num_blocks Number of blocks in the list
*/
void normalize_FFTW(float *blocks, int w, int num_blocks)
{
  const float ONE_DIV_2w = 1.0 / (2.0 * w);
  const float ONE_DIV_SQRT_2 = 1.0 / sqrtf(2);

  // Divide all coefficients by 2*w
  //#pragma omp parallel for shared(blocks)
  for (int i = 0; i < num_blocks * w * w; i++)
    blocks[i] *= ONE_DIV_2w;

#ifdef _OPENMP
#pragma omp parallel for shared(blocks) schedule(static)
#endif
  for (int b = 0; b < num_blocks; b++)
  {
    // {(i, j)} with i == 0 or j == 0
    for (int j = 0; j < w; j++)
    {
      int i = 0;
      blocks[b * w * w + j * w + i] *= ONE_DIV_SQRT_2;
    }
    for (int i = 0; i < w; i++)
    {
      int j = 0;
      blocks[b * w * w + j * w + i] *= ONE_DIV_SQRT_2;
    }
  }
}

/**
 * @brief Build a mask for invalid pixel. If mask(i, j) = true, the pixels will not be used.
 *
 * @param i_im : noisy image;
 * @param o_mask : will contain the mask for all pixel in the image size;
 * @param p_imSize : size of the image;
 * @param p_sizePatch : size of a patch.
 *
 * @return number of valid blocks.
 *
 **/
unsigned buildMask(CImage &i_im, int *o_mask,
                   unsigned Nx, unsigned Ny, unsigned w,
                   unsigned num_channels)
{
  unsigned count = 0;

  for (unsigned ij = 0; ij < Nx * Ny; ij++)
  {
    const unsigned j = ij / Nx;
    const unsigned i = ij - j * Nx;

    //! Look if the pixel is not to close to the boundaries of the image
    if (i < Nx - w + 1 && j < Ny - w + 1)
    {
      for (unsigned c = 0; c < num_channels; c++)
      {
        float *u = i_im.get_channel(c);

        //! Look if the square 2x2 of pixels is constant
        int invalid_pixel = (c == 0 ? 1 : o_mask[ij]);

        // Try to validate pixel
        if (fabs(u[ij] - u[ij + 1]) > 0.001f)
        {
          invalid_pixel = 0;
        }
        else if (fabs(u[ij + 1] - u[ij + Nx]) > 0.001f)
        {
          invalid_pixel = 0;
        }
        else if (fabs(u[ij + Nx] - u[ij + Nx + 1]) > 0.001f)
        {
          invalid_pixel = 0;
        }
        o_mask[ij] = invalid_pixel;
      }
    }
    else
    {
      o_mask[ij] = 1; // Not valid
    }

    if (o_mask[ij] == 0)
      count++;
  }

  return count;
}


/**
 * @brief Build a mask for saturated (invalid) pixel. If mask(i, j) = true, the pixels will not be used.
 *
 * @param i_im : noisy image;
 * @param o_mask : will contain the mask for all pixel in the image size;
 *
 * @return number of valid blocks.
 *
 **/
unsigned buildSaturatedMask(CImage &i_im, int *o_mask,
                    unsigned Nx, unsigned Ny, unsigned w,
                    unsigned num_channels) {
  unsigned count  = 0;
  
  //! Get the maximum value of the channel
  // float max_val = 0;
  std::vector<float> max_vals(num_channels);

  
  for (unsigned c = 0; c < num_channels; c++) {
    for (unsigned ij = 0; ij < Nx*Ny; ij++) {
      float *u = i_im.get_channel(c);
      max_vals[c] = max_vals[c] < u[ij] ? u[ij] : max_vals[c];
    }
  }

  memset(o_mask, 0, Nx*Ny*sizeof(*o_mask));

  for (unsigned j = 0; j < Ny; ++j) {
    for (unsigned i = 0; i < Nx; ++i) {
      const unsigned ij = j * Nx + i;

      for (unsigned c = 0; c < num_channels; c++) {
        float *u = i_im.get_channel(c);
        int invalid_pixel = (c == 0 ? 0 : o_mask[ij]);
        //! Look if the pixel is saturated
        o_mask[ij] = u[ij] >= max_vals[c] ? 1 : invalid_pixel;
      }
    }
  }

  for (unsigned j = 0; j < Ny; ++j) {
    for (unsigned i = 0; i < Nx; ++i) {
      const unsigned ij = j * Nx + i;

      if (i < Nx - w + 1 && j < Ny - w + 1) {
        unsigned accum = 0;
        for (size_t ii = 0; ii < w; ++ii) {
          for (size_t jj = 0; jj < w; ++jj) {
            accum += o_mask[ij + ii*Nx + jj];
          }
        }
        
        if (accum > 0) {
          o_mask[ij] = 1;
        } else {
          ++count;
        }
      } else {
        o_mask[ij] = 1;
      }
    }
  }


  // for (unsigned ij = 0; ij < Nx*Ny; ij++) {
  //   const unsigned j = ij / Nx;
  //   const unsigned i = ij - j * Nx;
    
  //   if (i < Nx - w + 1 && j < Ny - w + 1) {
  //     unsigned accum = 0;
  //     for (int ii = 0; ii < w; ++ii) {
  //       for (int jj = 0; jj < w; ++jj) {
  //         accum += o_mask[ij + ii*Nx + jj];
  //       }
  //     }
      
  //     if (accum > 0) {
  //       o_mask[ij] = 1;
  //     } else {
  //       ++count;
  //     }
  //   } else {
  //     o_mask[ij] = 1;
  //   }
  // }

  return count;
}

/**
 * @brief Build a mask of all pixels excluding right and bottom borders. If mask(i, j) = true, the pixels will not be used.
 *
 * @param o_mask : will contain the mask for all pixel in the image size;
 *
 * @return number of valid blocks.
 *
 **/
unsigned buildFullMask(int *o_mask,
                            unsigned Nx, unsigned Ny, unsigned w)
{

  for (size_t ij = 0; ij < Nx*Ny; ++ij) {
    const size_t j = ij / Nx;
    const size_t i = ij - j * Nx;
    o_mask[ij] = (i < Nx - w + 1 && j < Ny - w + 1) ? 0 : 1;
  }



  return (Nx - w + 1) * (Ny - w + 1);
}

/**
 * @brief Get the names of all .png or .bmp files in the given directory
 *
 * @param path : directory path
 *
 * @return an array of file names
 *
 **/
std::vector<std::string> get_filenames(const std::string &path)
{
  std::vector<std::string> files;
  struct dirent *entry;
  DIR *dir = opendir(path.c_str());

  if (dir == NULL)
  {
    return files;
  }

  while ((entry = readdir(dir)) != NULL)
  {
    std::string fname(entry->d_name);
    if (fname.size() > 4)
    {
      if (fname.substr(fname.size() - 4, fname.size()).compare(".png") == 0 ||
          fname.substr(fname.size() - 4, fname.size()).compare(".bmp") == 0)
      {
        files.push_back(path + "/" + fname);
      }
    }
  }
  closedir(dir);

  return files;
}

//! Computes the set of variances of all bins, computed from the high-frequency coefficients of the blocks
// from multiple images.
/*!
  \param *VH Output sets of variances of all bins, of size num_bins*w*w
  \param block_ids_selected List of pointers to the blocks
  \param inputs List of pointers to the blocks
  \param w Block side
  \param T Number of low-freq coefficients, excluding DC
  \param ch Channel of the computed variances
  \return Length of the returned variances list
*/
int compute_VH(double *VH, const std::vector<BlockId> &block_ids_selected,
               const std::vector<std::string> &inputs, int w, int T, int ch)
{
  std::vector<std::vector<BlockId>> block_ids_by_image(inputs.size());
  for (const BlockId block_id : block_ids_selected)
  {
    block_ids_by_image[block_id.image].push_back(block_id);
  }

  int VH_count = 0;
  for (uint32_t image_idx = 0; image_idx < inputs.size(); ++image_idx)
  {
    PRINT_VERBOSE("Loading again %s for computing VH\n", inputs[image_idx].c_str());

    // load image
    CImage input;
    input.load(inputs[image_idx].c_str());

    int Nx = input.get_width();
    // int Ny = input.get_height();

    // Initialize blocks
    const std::vector<BlockId> &block_ids = block_ids_by_image[image_idx];
    uint32_t num_blocks = block_ids.size();

    PRINT_VERBOSE("num_blocks: %d\n", num_blocks);

    float *blocks = new float[num_blocks * w * w];

    int nbTable[2] = {w, w};
    int nembed[2] = {w, w};

#ifdef _OPENMP
    fftwf_plan_with_nthreads(omp_get_num_procs());
#endif

    fftwf_r2r_kind kindTable[2] = {FFTW_REDFT10, FFTW_REDFT10};

    fftwf_plan fft_plan = fftwf_plan_many_r2r(2, nbTable, num_blocks, blocks,
                                              nembed, 1, w * w, blocks, nembed,
                                              1, w * w, kindTable, FFTW_ESTIMATE);

    const uint32_t w2 = w * w;
    float *image = input.get_channel(ch);


#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t q = 0; q < num_blocks; ++q)
    {
      for (int j = 0; j < w; ++j)
      {
        for (int i = 0; i < w; ++i)
        {
          int block_origin = block_ids[q].block_origin;
          blocks[q * w2 + j * w + i] = image[block_origin + j * Nx + i];
        }
      }
    }

    // Compute 2D-DCT of all the blocks
    //
    // Transform blocks with FFTW
    fftwf_execute_r2r(fft_plan, blocks, blocks);

    // Normalize FFTW output
    normalize_FFTW(blocks, w, num_blocks);

    // store VH for each bin
    VH_count = 0;
    for (int q = 0; q < w * w; ++q)
    {
      int j = q / w;
      int i = q - j * w;
      if (i + j >= T)
      {
        for (uint32_t k = 0; k < num_blocks; ++k)
        {
          const BlockId &block_id = block_ids[k];
          VH[block_id.bin * w2 + VH_count] += pow(blocks[k * w2 + q], 2);
        }
        ++VH_count;
      }
    }

    delete[] blocks;

    fftwf_destroy_plan(fft_plan);
  }

  return VH_count;
}

//! Ponomarenko et al. AVIRIS noise estimation algorithm.
/*!
  \param argc Number of arguments of the program
  \param **argv Arguments of the program
*/
void algorithm(int argc, char **argv)
{
  vector<OptStruct *> options;
  vector<ParStruct *> parameters;
  //
  OptStruct owin = {"w:", 8, "8", NULL, "Block side"};
  options.push_back(&owin);
  OptStruct opercentile = {"p:", 1, "0.005", NULL, "Percentile"};
  options.push_back(&opercentile);
  // OptStruct ore = {"r", 0, NULL, NULL, "Flag to remove equal pixels"};
  // options.push_back(&ore);
  OptStruct obins = {"b:", 0, "0", NULL, "Number of bins"};
  options.push_back(&obins);
  OptStruct oD = {"D:", 7, "7", NULL, "Filtering distance"};
  options.push_back(&oD);
  OptStruct ofiltercurve = {"g:", 5, "5", NULL, "Filter curve iterations"};
  options.push_back(&ofiltercurve);
  OptStruct omeanMethod = {"m:", 2, "2", NULL, "Mean computation method"};
  options.push_back(&omeanMethod);
  OptStruct oremoveSaturate = {"s", 0, NULL, NULL, "Flag to remove saturated pixels"};
  options.push_back(&oremoveSaturate);  
  OptStruct ochannel = {"c:", 3, "3", NULL, "Number of channels"};
  options.push_back(&ochannel);

  ParStruct pinput = {"input", NULL, "input image directory"};
  parameters.push_back(&pinput);
  //
  if (!parsecmdline("ponomarenko", "Ponomarenko SD noise estimation algorithm",
                    argc, argv, options, parameters))
  {
    printf("\n");
    printf("(c) 2012 Miguel Colom. Under license GNU GPL.\n");
    printf("http://mcolom.perso.math.cnrs.fr/\n");
    printf("\n");
    exit(-1);
  }

  // Read parameters
  int w = atoi(owin.value);
  int T = get_T(w);
  float p = atof(opercentile.value);
  int num_bins = atoi(obins.value);
  int D = atoi(oD.value);
  int curve_filter_iterations = atoi(ofiltercurve.value);
  int mean_method = atoi(omeanMethod.value);
  // bool remove_equal_pixels_blocks = ore.flag;
  bool remove_saturated_pixels_blocks = oremoveSaturate.flag;
  int num_channels = atoi(ochannel.value);

// Parallelization config
#ifdef _OPENMP
  omp_set_num_threads(omp_get_num_procs());
#endif

  CFramework::set_verbose(false);

  // Custom percentile or given by the user?
  // bool custom_percentile = is_custom_percentile(opercentile.value);

  // Load input image
  std::string path(pinput.value);
  std::vector<std::string> inputs = get_filenames(path);

  // TODO: prepare statistics variable
  std::vector<std::vector<float>> means_all(num_channels);
  std::vector<std::vector<float>> VL_all(num_channels);
  std::vector<std::vector<BlockId>> block_ids_all(num_channels);

  // float *vmeans = NULL;
  // float *vstds = NULL;
  std::vector<float> vmeans;
  std::vector<float> vstds;

  for (uint32_t input_idx = 0; input_idx < inputs.size(); ++input_idx)
  {

    // TODO: input multiple images
    PRINT_VERBOSE("Processing %s \n", inputs[input_idx].c_str());

    CImage input;
    input.load(inputs[input_idx].c_str());

    // Get image properties
    int Nx = input.get_width();
    int Ny = input.get_height();
    int input_channels = input.get_num_channels();

    if (input_channels != num_channels)
    {
      printf("%s has a channel number (%d) incompatible with the required channel number (%d). \n"
             "Skipped.\n",
             inputs[input_idx].c_str(),
             input_channels, num_channels);

      continue;
    }

    int total_blocks = (Nx - w + 1) * (Ny - w + 1); // Number of overlapping blocks

    // Create equal pixels mask
    int *mask_all = NULL;
    int num_blocks;


    // int *mask_saturated;
    mask_all = new int[Nx * Ny];
    if (remove_saturated_pixels_blocks) {
      num_blocks = buildSaturatedMask(input, mask_all, Nx, Ny, w, num_channels);
    } 
    // else if (remove_equal_pixels_blocks)
    // {
    //   mask_all = new int[Nx * Ny];
    //   num_blocks = buildMask(input, mask_all, Nx, Ny, w, num_channels);
    // }
    else
    {
      num_blocks = buildFullMask(mask_all, Nx, Ny, w);
    }


    if (input_idx == 0)
    {
      // reserve space for global storage
      for (int ch = 0; ch < num_channels; ++ch)
      {
        means_all[ch].reserve(inputs.size() * total_blocks);
        VL_all[ch].reserve(inputs.size() * total_blocks);
        block_ids_all[ch].reserve(inputs.size() * total_blocks);
      }

      // Set number of bins
      if (num_bins <= 0)
        num_bins = Nx * Ny / 42000;
      if (num_bins <= 0)
        num_bins = 1; // Force at least one bin

      // Initialize the arrays for final means and noise estimations
      // vmeans = new float[num_channels * num_bins];
      // vstds = new float[num_channels * num_bins];
      vmeans.resize(num_channels * num_bins);
      vstds.resize(num_channels * num_bins);
    }

    // Compute delta and theta
    CFramework *fw = CFramework::get_framework();
    float *delta = fw->create_array(w * w);
    int theta = compute_delta(delta, w, T);

    std::vector<float> means(num_blocks);

    std::vector<BlockId> block_ids(num_blocks);
    float *blocks = new float[num_blocks * w * w];

    // Init FFTW threads
    fftwf_init_threads();

    int nbTable[2] = {w, w};
    int nembed[2] = {w, w};

#ifdef _OPENMP
    fftwf_plan_with_nthreads(omp_get_num_procs());
#endif

    fftwf_r2r_kind kindTable[2] = {FFTW_REDFT10, FFTW_REDFT10};

    fftwf_plan fft_plan = fftwf_plan_many_r2r(2, nbTable, num_blocks, blocks,
                                              nembed, 1, w * w, blocks, nembed,
                                              1, w * w, kindTable, FFTW_ESTIMATE);
    
    unsigned* valid_coords = new unsigned[num_blocks];
    make_valid_coords(mask_all, valid_coords, Nx, Ny);

    // Process each channel
    for (int ch = 0; ch < num_channels; ch++)
    {
      float *u = input.get_channel(ch);
      
      
      read_all_valid_blocks(blocks, u, Nx, Ny, w, num_blocks, mask_all);

      // Compute means
      compute_means(means.data(), blocks, w, num_blocks);

      // Compute 2D-DCT of all the blocks
      //
      // Transform blocks with FFTW
      fftwf_execute_r2r(fft_plan, blocks, blocks);

      // Normalize FFTW output
      normalize_FFTW(blocks, w, num_blocks);

      // Create a list of pointers of the groups
      float **blocks_ptr = new float *[num_blocks];
      for (int i = 0; i < num_blocks; i++)
        blocks_ptr[i] = &blocks[i * w * w];

      // Create the global indices for blocks
      for (int i = 0; i < num_blocks; ++i)
      {
        block_ids[i].image = input_idx;
        // block_ids[i].block_origin = i;
        block_ids[i].block_origin = valid_coords[i];
        block_ids[i].bin = -1;
      }
        // cout << "valid_coords[" << i << "]" << endl;

      // Compute VL
      std::vector<float> VL(num_blocks);
      compute_VL(VL.data(), num_blocks, w, delta, blocks_ptr, theta);

      // store the means, VL and block indices
      VL_all[ch].insert(VL_all[ch].end(), VL.begin(), VL.end());
      means_all[ch].insert(means_all[ch].end(), means.begin(), means.end());
      block_ids_all[ch].insert(block_ids_all[ch].end(), block_ids.begin(), block_ids.end());

      /*
      // Create histogram according to the means
      CHistogram<float*> histo = CHistogram<float*>(num_bins,
                                                    blocks_ptr,
                                                    means.data(),
                                                    num_blocks);

      // Process each bin
      #ifdef _OPENMP
      #pragma omp parallel for shared(vmeans, vstds, histo) schedule(static)
      #endif
      for (int bin = 0; bin < num_bins; bin++) {
        int elems_bin = histo.get_num_elements_bin(bin);

        float **block_ptr_bin = histo.get_data_bin(bin);

        float *VL = new float[elems_bin];

        // Compute VL
        compute_VL(VL, elems_bin, w, delta, block_ptr_bin, theta);

        
      }*/

      delete[] blocks_ptr;
    }
    if (mask_all != NULL)
      delete[] mask_all;

    delete[] blocks;
    delete[] valid_coords;
  } // end for each input

  for (int ch = 0; ch < num_channels; ++ch)
  {
    assert(VL_all[ch].size() == means_all[ch].size());
    assert(block_ids_all[ch].size() == means_all[ch].size());

    // Create histogram according to the means
    CHistogram<BlockId> histo_block_id = CHistogram<BlockId>(num_bins,
                                                             block_ids_all[ch].data(),
                                                             means_all[ch].data(),
                                                             means_all[ch].size());

    // Create another histogram according to the means to store VL
    CHistogram<float> histo_VL = CHistogram<float>(num_bins,
                                                   VL_all[ch].data(),
                                                   means_all[ch].data(),
                                                   means_all[ch].size());

    // print_arr(means_all[ch].data(), means_all[ch].size());

    // store the blocks used for VH computation
    // The relationship is: image <-> block <-> bin
    std::vector<BlockId> block_ids_selected;

    // for (int i = 0; i < block_ids_all[ch].size(); ++i) {
    //   auto e = block_ids_all[ch][i];
    //   if (e.image == 1) continue;
    //   printf("image, block, bin = %d, %d, %d\n", e.image, e.block_origin, e.bin);

    //   printf("mean, VL = %f, %f\n", means_all[ch][i], VL_all[ch][i]);

    // }

    // reserve more space for storage
    block_ids_selected.reserve(block_ids_all[ch].size() * p * 1.1);

    for (int bin = 0; bin < num_bins; ++bin)
    {
      int elems_bin = histo_block_id.get_num_elements_bin(bin);

      PRINT_VERBOSE("elems_bin: %d\n", elems_bin);

      BlockId *block_ids = histo_block_id.get_data_bin(bin);
      float *VL = histo_VL.get_data_bin(bin);

      // Use directly the fixed quantile
      int K = elems_bin * p;

      // if (custom_percentile)
      //   K = elems_bin * p;
      // else // Using Ponomarenko's article loop
      //   K = get_optimal_K_ponom_orig(elems_bin, VL);

      int *indices_VL = new int[elems_bin];
      argsort(VL, indices_VL, elems_bin);

      // for (int i = 0; i < K; ++i) {
      //   printf("VL: %f, image: %f", VL[indices_VL[i]], )
      // }

      for (int i = 0; i < K; ++i)
      {
        BlockId block_id = block_ids[indices_VL[i]];
        block_id.bin = bin;
        block_ids_selected.push_back(block_id);
      }

      float *means_bin = histo_block_id.get_datal_bin(bin);
      float bin_mean = get_top_K_mean(mean_method, K, indices_VL, means_bin);

      vmeans[ch * num_bins + bin] = bin_mean;

      delete[] indices_VL;
    }

    std::vector<double> VH(num_bins * w * w, 0);
    int VH_count = compute_VH(VH.data(), block_ids_selected, inputs, w, T, ch);

    PRINT_VERBOSE("VH_count: %d\n", VH_count);

    for (int bin = 0; bin < num_bins; ++bin)
    {
      int K = histo_block_id.get_num_elements_bin(bin) * p;
      // Normalize the VH for each bin
      for (int i = 0; i < VH_count; ++i)
      {
        VH[bin * w * w + i] /= K;
      }

      PRINT_VERBOSE("VH of bin %d: \n ", bin);
      for (int c = 0; c < VH_count; ++c)
      {
        PRINT_VERBOSE("%f, ", VH[bin * w * w + c]);
      }
      PRINT_VERBOSE("\n");

      float tilde_sigma = sqrt(median(&VH[bin * w * w], VH_count));
      vstds[ch * num_bins + bin] = tilde_sigma;

      PRINT_VERBOSE("K: %d\n", int(histo_block_id.get_num_elements_bin(bin) * p));
    }

  }

  PRINT_VERBOSE("vstds:\n");
  for (int bin = 0; bin < num_bins; ++bin)
  {
    for (int ch = 0; ch < num_channels; ++ch)
    {
      PRINT_VERBOSE("%f ", vstds[ch * num_bins + bin]);
    }
    PRINT_VERBOSE("\n");
  }

  // Filter noise curve
  float *new_std_control = new float[num_bins * num_channels];
  copy(new_std_control, vstds.data(), num_channels * num_bins);
  //
  for (int ch = 0; ch < num_channels; ch++)
    for (int filt_iter = 0; filt_iter < curve_filter_iterations; filt_iter++)
    {
      bool allow_up = (filt_iter < 3);
      filter_curve(&vmeans[ch * num_bins], &new_std_control[ch * num_bins],
                   num_bins,
                   &new_std_control[ch * num_bins],
                   D, allow_up);
    }

  // Print results
  for (int bin = 0; bin < num_bins; bin++)
  {
    // Means
    for (int ch = 0; ch < num_channels; ch++)
      printf("%f  ", vmeans[ch * num_bins + bin]);

    // Standard deviations
    for (int ch = 0; ch < num_channels; ch++)
      printf("%f  ", new_std_control[ch * num_bins + bin]);
    //
    printf("\n");
  }

  // FFTW Cleanup

  fftwf_cleanup_threads();
  fftwf_cleanup();

  delete[] new_std_control;
  // if (vmeans != NULL)
  //   delete[] vmeans;
  // if (vstds != NULL)
  //   delete[] vstds;
}
