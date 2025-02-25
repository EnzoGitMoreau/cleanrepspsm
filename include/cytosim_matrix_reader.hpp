#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "matrix44.h"
#include "sparmatsymblk.h"
#ifndef BLOCKSIZE
	#define BLOCKSIZE 4
#endif 
#ifdef RSB
#include <rsb.hpp>
#endif
#ifdef MKL
    #include <mkl.h>
    #include <mkl_types.h>
    #include <mkl_spblas.h>
    #endif

class CytMatrixReader
{
private:
	std::vector<Matrix44*> blocks;
	short** isCreatedBlock;
	short** indexInVector;
	int matrixSize;
	int internMatSize;
	int numberofEntries;
	int number_of_block_per_line;
public:
	std::vector<std::string> splitStringByWhitespace(const std::string& str) {
	std::vector<std::string> result;
	std::istringstream iss(str);
	std::string token;
	
	

	while (iss >> token) 
	{
			result.push_back(token);
	}

		return result;
	}
	
	#ifndef RSB
	#ifndef MKL
		CytMatrixReader(std::string matrix_name,std::string vector_name, SparMatSymBlk* matrix, real** vector,int nb_threads)
		{
			//Reading Matrix
			std::ifstream inputFile(matrix_name);
			if(!inputFile)
			{
				VLOGe("Could not read "+matrix_name);	
			}
			else
			{
				
				std::string ligne;
				int i =0;
				std::istringstream iss;
				
				std::vector<std::string> splitLine;
				
				while (getline(inputFile, ligne)) 
				{ 
					splitLine = splitStringByWhitespace(ligne);
					if(i==1)
					{
						matrixSize = std::stoi(splitLine[0]);
						numberofEntries = std::stoi(splitLine[2]);
						VLOG("Cytosim's Matrix informations:");
						VLOG("Matrix size: "+CONV(matrixSize)+" Number of entries: "+CONV(numberofEntries));
						VLOG("Reading and constructing Matrix");
						if(matrixSize%(BLOCKSIZE*nb_threads) !=0)
						{
						matrix->resize((int(matrixSize/(BLOCKSIZE*nb_threads))+1)*(BLOCKSIZE*nb_threads));
						internMatSize = (int(matrixSize/(BLOCKSIZE*nb_threads))+1)*(BLOCKSIZE*nb_threads);
						}
						else
						{
							matrix->resize(matrixSize);
							internMatSize = matrixSize;
						}
					}
					if(i>1)
					{
				
						std::istringstream iss(splitLine[2]);
						
						int i = std::stoi(splitLine[0]);
					
						int j = std::stoi(splitLine[1]);
					
						real value;
						iss>>value;
						real& value_in_matrix = matrix->element(i,j);
						value_in_matrix = (int) value;
					
						
					}
					i++;
				}
			}

			std::ifstream inputFileV(vector_name);
			if(!inputFileV)
			{
				VLOGe("Could not read "+vector_name);	
			}
			else
			{
				
				std::string ligne;
				int i =0;
				std::istringstream iss;
				
				std::vector<std::string> splitLine;
				
				while (getline(inputFileV, ligne)) 
				{ 
					splitLine = splitStringByWhitespace(ligne);
					if(i==1)
					{
						int vector_size = std::stoi(splitLine[0]);
						if(vector_size != matrixSize)
						{
							numberofEntries = 0;
							matrixSize = 0;
							VLOGe("Error: Matrix size and vector size are different");
						}
						else
						{
						
						numberofEntries = vector_size;
						VLOG("Cytosim's Vector informations:");
						VLOG("Vector size: "+CONV(vector_size)+" Number of entries: "+CONV(numberofEntries));
						VLOG("Reading and constructing vector");
						*vector = (real*) malloc(sizeof(real)*internMatSize);//Maybe
						}
					}
					if(i>1)
					{
						std::istringstream iss(splitLine[1]);
						int ind = std::stoi(splitLine[0]);
						real val;
						iss>> val;
						val = (int) val;
						(*vector)[ind] = val;
					
						
					}
					i++;
				}
			}
			//Reading vector

			

		}
	#endif 
		#ifdef MKL
		CytMatrixReader(std::string matrix_name, std::string vector_name, SparMatSymBlk* matrix,sparse_matrix_t** mtx_ptr, real** vector,int nb_threads)
	{
	
		
		std::ifstream inputFile(matrix_name);
		if(!inputFile)
		{
			VLOGe("Could not read "+matrix_name);	
		}
		else
		{
			
			std::string ligne;
			int i =0;
			std::istringstream iss;
			short** already_put_index;
			std::vector<std::string> splitLine;
			//Reading matrix
			while (getline(inputFile, ligne)) 
			{ 
				splitLine = splitStringByWhitespace(ligne);
				if(i==1)
				{
					matrixSize = std::stoi(splitLine[0]);
					numberofEntries = std::stoi(splitLine[2]);
					VLOG("Cytosim's Matrix informations:");
					VLOG("Matrix size: "+CONV(matrixSize)+" Number of entries: "+CONV(numberofEntries));
					VLOG("Reading and constructing Matrix");
					if(matrixSize%(BLOCKSIZE*nb_threads) !=0)
					{
					matrix->resize((int(matrixSize/(BLOCKSIZE*nb_threads))+1)*(BLOCKSIZE*nb_threads));
					internMatSize = (int(matrixSize/(BLOCKSIZE*nb_threads))+1)*(BLOCKSIZE*nb_threads);
					}
					else
					{
						matrix->resize(matrixSize);
						internMatSize = matrixSize;
					}
					
					//const rsb_coo_idx_t nrA { matrixSize }, ncA { matrixSize };
					//*mtx_ptr = new rsb::RsbMatrix<double>(nrA, ncA);
					
					
				}
				if(i>1)
				{
					std::istringstream iss(splitLine[2]);
					int i = std::stoi(splitLine[0]);
					int j = std::stoi(splitLine[1]);
					real value;
					iss>>value;
					real& value_in_matrix = matrix->element(i,j);
					value =  value;
					value_in_matrix =value;//temp
					
					
					
					

					
				
					
				}
				i++;
			}
			
			VLOG("Done reading input Matrix");
			//Reading vector
			std::ifstream inputFileV(vector_name);
			if(!inputFileV)
			{
				VLOGe("Could not read "+vector_name);
			}
			else
			{
				
				std::string ligne;
				int i =0;
				std::istringstream iss;
				
				std::vector<std::string> splitLine;
				
				while (getline(inputFileV, ligne)) 
				{ 
					
					
					if(i==1)
					{
						int vector_size = std::stoi(ligne);
						if(vector_size != matrixSize)
						{
							numberofEntries = 0;
							matrixSize = 0;
							VLOGe("Error: Matrix size and vector size are different");
						}
						else
						{
						
						numberofEntries = vector_size;
						VLOG("Cytosim's Vector informations:");
						VLOG("Vector size: "+CONV(vector_size)+" Number of entries: "+CONV(numberofEntries));
						VLOG("Reading and constructing vector");
						*vector = (real*) malloc(sizeof(real)*internMatSize);//Maybe
						}
					}
					if(i>1)
					{
						std::istringstream iss(ligne);
						real value;
						iss>>value;
						value = (int) value;
						(*vector)[i-2] = value;
					
						
					}
					i++;
				}
			}
			}
	
		

	}
		#endif
	#endif
	#ifdef RSB
	CytMatrixReader(std::string matrix_name, std::string vector_name, SparMatSymBlk* matrix,rsb::RsbMatrix<double>** mtx_ptr, real** vector,int nb_threads)
	{
		rsb::RsbLib rsblib;
		
		std::ifstream inputFile(matrix_name);
		if(!inputFile)
		{
			VLOGe("Could not read "+matrix_name);	
		}
		else
		{
			
			std::string ligne;
			int i =0;
			std::istringstream iss;
			short** already_put_index;
			std::vector<std::string> splitLine;
			//Reading matrix
			while (getline(inputFile, ligne)) 
			{ 
				splitLine = splitStringByWhitespace(ligne);
				if(i==1)
				{
					matrixSize = std::stoi(splitLine[0]);
					numberofEntries = std::stoi(splitLine[2]);
					VLOG("Cytosim's Matrix informations:");
					VLOG("Matrix size: "+CONV(matrixSize)+" Number of entries: "+CONV(numberofEntries));
					VLOG("Reading and constructing Matrix");
					if(matrixSize%(BLOCKSIZE*nb_threads) !=0)
					{
					matrix->resize((int(matrixSize/(BLOCKSIZE*nb_threads))+1)*(BLOCKSIZE*nb_threads));
					internMatSize = (int(matrixSize/(BLOCKSIZE*nb_threads))+1)*(BLOCKSIZE*nb_threads);
					}
					else
					{
						matrix->resize(matrixSize);
						internMatSize = matrixSize;
					}
					
					const rsb_coo_idx_t nrA { matrixSize }, ncA { matrixSize };
					*mtx_ptr = new rsb::RsbMatrix<double>(nrA, ncA);
					
					
				}
				if(i>1)
				{
					std::istringstream iss(splitLine[2]);
					int i = std::stoi(splitLine[0]);
					int j = std::stoi(splitLine[1]);
					real value;
					iss>>value;
					real& value_in_matrix = matrix->element(i,j);
					value =  value;
					value_in_matrix =value;//temp
					double errval = (*mtx_ptr)->get_val(i,j);
					if(!errval)
					{(*mtx_ptr)->set_val(value, i,j);
						if(i!=j)
						{
						(*mtx_ptr)->set_val(value, j,i);
						}}
					else
					{
						std::cout<<"error";
					}
					
					
					

					
				
					
				}
				i++;
			}
			(*mtx_ptr)->close();
			VLOG("Done reading input Matrix");
			//Reading vector
			std::ifstream inputFileV(vector_name);
			if(!inputFileV)
			{
				VLOGe("Could not read "+vector_name);
			}
			else
			{
				
				std::string ligne;
				int i =0;
				std::istringstream iss;
				
				std::vector<std::string> splitLine;
				
				while (getline(inputFileV, ligne)) 
				{ 
					
					
					if(i==1)
					{
						int vector_size = std::stoi(ligne);
						if(vector_size != matrixSize)
						{
							numberofEntries = 0;
							matrixSize = 0;
							VLOGe("Error: Matrix size and vector size are different");
						}
						else
						{
						
						numberofEntries = vector_size;
						VLOG("Cytosim's Vector informations:");
						VLOG("Vector size: "+CONV(vector_size)+" Number of entries: "+CONV(numberofEntries));
						VLOG("Reading and constructing vector");
						*vector = (real*) malloc(sizeof(real)*internMatSize);//Maybe
						}
					}
					if(i>1)
					{
						std::istringstream iss(ligne);
						real value;
						iss>>value;
						value = (int) value;
						(*vector)[i-2] = value;
					
						
					}
					i++;
				}
			}
			}
	
		

	}
	#endif
	
	int problemSize()
	{
		return matrixSize;
	}
	
	
};