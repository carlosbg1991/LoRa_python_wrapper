/* file : lora.i */
  
/* name of module to use*/
%module lora
%{ 
    /* Every thing in this file is being copied in  
     wrapper file. We include the C header file necessary 
     to compile the interface */
    #include "decoder.h"
%} 
  
/* explicitly list functions and variables to be interfaced
   or if we want to interface all functions then we can simply 
   include header file like this - %include "decoder.h" */
%include "decoder.h"
