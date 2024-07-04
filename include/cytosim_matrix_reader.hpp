#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "matrix44.h"
#include "sparmatsymblk.h"

#ifdef RSB
#include <rsb.hpp>
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
	CytMatrixReader(std::string matrix_name,std::string vector_name, SparMatSymBlk* matrix, real** vector,int nb_threads)
	{
		//Reading Matrix
		std::ifstream inputFile(matrix_name);
		if(!inputFile)
		{
			std::cout<<"\nCouldn't read "<<matrix_name<<"\n";	
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
					std::cout<<"\n[INFO] Cytosim's Matrix informations: \n";
					std::cout << "[INFO] Matrix size: "<<matrixSize <<" Number of entries: "<<numberofEntries;
					std::cout <<"\n[INFO] Reading and constructing Matrix";
					if(matrixSize%(4*nb_threads) !=0)
					{
					matrix->resize((int(matrixSize/(4*nb_threads))+1)*(4*nb_threads));
					internMatSize = (int(matrixSize/(4*nb_threads))+1)*(4*nb_threads);
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
			std::cout<<"\nCouldn't read "<<vector_name<<"\n";	
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
						std::cout<<"Error: Matrix size and vector size are different";
					}
					else
					{
					
					numberofEntries = vector_size;
					std::cout<<"\n[INFO] Cytosim's Vector informations: \n";
					std::cout << "[INFO] Vector size: "<<vector_size<<" Number of entries: "<<numberofEntries;
					std::cout <<"\n[INFO] Reading and constructing vector";
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
	#ifdef RSB
	CytMatrixReader(std::string matrix_name, std::string vector_name, SparMatSymBlk* matrix,rsb::RsbMatrix<double>** mtx_ptr, real** vector,int nb_threads)
	{
		rsb::RsbLib rsblib;
		
		std::ifstream inputFile(matrix_name);
		if(!inputFile)
		{
			std::cout<<"\nCouldn't read "<<matrix_name<<"\n";	
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
					std::cout<<"\n[INFO] Matrix Market's Matrix informations: \n";
					std::cout << "[INFO] Matrix size: "<<matrixSize <<" Number of entries: "<<numberofEntries;
					std::cout <<"\n[INFO] Reading and constructing Matrix";
					if(matrixSize%(4*nb_threads) !=0)
					{
					matrix->resize((int(matrixSize/(4*nb_threads))+1)*(4*nb_threads));
					internMatSize = (int(matrixSize/(4*nb_threads))+1)*(4*nb_threads);
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
					value_in_matrix = value;
					double errval = (*mtx_ptr)->get_val(i,j);
					if(!errval)
					{(*mtx_ptr)->set_val(value, i,j);
						if(i!=j)
						{
						(*mtx_ptr)->set_val(value, j,i);
						}}
					
					
					

					
				
					
				}
				i++;
			}
			(*mtx_ptr)->close();
			std::cout<<"[INFO] Done reading input Matrix \n";
			//Reading vector
			std::ifstream inputFileV(vector_name);
			if(!inputFileV)
			{
				std::cout<<"\nCouldn't read "<<vector_name<<"\n";	
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
							std::cout<<"Error: Matrix size and vector size are different";
						}
						else
						{
						
						numberofEntries = vector_size;
						std::cout<<"\n[INFO] Cytosim's Vector informations: \n";
						std::cout << "[INFO] Vector size: "<<vector_size<<" Number of entries: "<<numberofEntries;
						std::cout <<"\n[INFO] Reading and constructing vector";
						*vector = (real*) malloc(sizeof(real)*internMatSize);//Maybe
						}
					}
					if(i>1)
					{
						std::istringstream iss(splitLine[1]);
						int ind = std::stoi(splitLine[0]);
			
						iss>>(*vector)[ind];
					
						
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