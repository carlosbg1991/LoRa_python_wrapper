/* file : lora.i */
/* -*- Python -*- #
*
* Copyright 2020 Carlos Bocanegra.
* The Genesys Lab, Northeastern University, Boston, MA
*
*/
  
/* name of module to use*/
%module lora_id
%include "std_string.i"
%include "std_complex.i"
%include "std_vector.i"
%include <std_list.i>
%include "cpointer.i"
%include "stdint.i"  // handles unisnged ints and more

// Instantiate templates used by module
%template(ComplexList) std::list<std::complex<float>>;
%template(ComplexVector) std::vector<std::complex<float>>;
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(FloatVector) std::vector<float>;

%inline %{
    /* Note: double[4][4] is equivalent to a pointer to an array double (*)[4] */
    double (*new_mat44())[4] {
        return (double (*)[4]) malloc(16*sizeof(double));
    }
    void free_mat44(double (*x)[4]) {
        free(x);
    }
    void mat44_set(double x[4][4], int i, int j, double v) {
        x[i][j] = v;
    }
    double mat44_get(double x[4][4], int i, int j) {
        return x[i][j];
    }
    std::complex<float>* new_grcomplex(int size) {
    	return (std::complex<float>*) malloc(size*sizeof(std::complex<float>));
    }
    void set_grcomplex(std::complex<float>* input, int index, std::complex<float> value){
    	input[index] = value;
    }
    void set_grcomplex_global(std::complex<float>* input, std::vector<std::complex<float>> vector_in){
    	input[0] = vector_in[0];
    }
    std::complex<float> get_grcomplex(std::complex<float>* input, int index){
    	return input[index];
    }
    unsigned int conv_unsigned_int(int input){
    	return (unsigned int) input;
    }
	void free_grcomplex(std::complex<float>* my_gr){
		free(my_gr);
	}
	int* create_int_ptr(int idx) {
    	int* my_ptr = (int*) malloc(sizeof(int));
    	my_ptr[0] = idx;
    	return my_ptr;
    }
    int get_int(int* my_int) {
    	return my_int[0];
    }
    struct loraphy_header_py {
	    uint8_t length;
	    uint8_t crc_msn : 4;
	    uint8_t has_mac_crc : 1;
	    uint8_t cr : 3;
	    uint8_t crc_lsn : 4;
	    uint8_t reserved : 4;
	};
	typedef loraphy_header_py d_phdr_py;
	struct DecoderState {
	    enum {
	        DETECT,
	        SYNC,
	        FIND_SFD,
	        PAUSE,
	        DECODE_HEADER,
	        DECODE_PAYLOAD,
	        STOP
	    };
	};
	DecoderState* create_decoderstate_ptr() {
    	DecoderState* d_state_ptr = (DecoderState*) malloc(sizeof(DecoderState));
    	return d_state_ptr;
    }
%}

%pythoncode %{
	d_phdr_py = loraphy_header_py
	d_state_py = DecoderState
%}

%{ 
    /* Every thing in this file is being copied in  
     wrapper file. We include the C header file necessary 
     to compile the interface */
    #include "decoder.h"
    #include "channelizer.h"
    #include "fractional_resampler.h"
%} 
  
/* explicitly list functions and variables to be interfaced
   or if we want to interface all functions then we can simply 
   include header file like this - %include "decoder.h" */
%include "channelizer.h"
%include "decoder.h"
%include "fractional_resampler.h"

