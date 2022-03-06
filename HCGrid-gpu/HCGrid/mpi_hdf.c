// An example to perform parallel reads from an existing HDF5 file.
// Command line inputs: filename, dataset name, dimX, dimY (2D dataset dimension)
// For this example, the number of processes should match the number of datasets read. 
 
#include "hdf5.h"
#include "stdlib.h"
#include <mpi.h>
#include <string.h>

#define DIMS   2

int
main (int argc, char **argv) {
    
    hid_t       fileId, datasetId1,datasetId2,datasetId3,datasetId4;
    char *v1;
    char *v2;
    char *v3;
    char *v4;
    hid_t       mspace,fspace;      
    hsize_t     dimsf[2];                 
    hsize_t     dimsm[2];
    int         *data, i;                   
    hid_t	plistId;             
    herr_t      status;

    int size, rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    // Initialize the MPI environment and get size and rank
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);  

    if ((argc < 5) || (argc > 8)) {
      printf("Usage error: please use mpirun -np <procs> ./<exe> <filename> <dataset_dimx> <dataset_dimy> <dataset name1> <..name2> <upto name4>\n");
      MPI_Finalize(); 
      return -1;
    }
    if (size > 4) {
      printf("Number of processes > 4 is not supported\n");
      MPI_Finalize();
      return -1;
    } 
    
    dimsf[0] = atoi(argv[2]);
    dimsf[1] = atoi(argv[3]);

    v1 = (char *)malloc(sizeof(char)*4);
    v2 = (char *)malloc(sizeof(char)*4);
    v3 = (char *)malloc(sizeof(char)*4);
    v4 = (char *)malloc(sizeof(char)*4);
   
    sprintf(v1,"%s",argv[4]);
    if (argc == 6) {
      sprintf(v2,"%s",argv[5]);
    } else if (argc == 7) {
      sprintf(v2,"%s",argv[5]);
      sprintf(v3,"%s",argv[6]);
    } else if (argc == 8) {
      sprintf(v2,"%s",argv[5]);
      sprintf(v3,"%s",argv[6]);
      sprintf(v4,"%s",argv[7]);
    }

    // Dimension of both memory space and filespace are the same for this example
    dimsm[0] = dimsf[0]; dimsm[1] = dimsf[1];

    // Allocate memory for the data being read
    data = (int *) malloc(sizeof(int)*dimsm[0]*dimsm[1]);

    //  Create file access property list to perform parallel I/O
    plistId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plistId, comm, info);

    // Open the file. Note that this is a collective operation
    fileId = H5Fopen(argv[1], H5F_ACC_RDONLY, plistId);
    H5Pclose(plistId);
   

    
    // Create the dataspace for the dataset.
    mspace = H5Screate_simple(DIMS, dimsf, NULL); 
    fspace = H5Screate_simple(DIMS, dimsm, NULL); 


    // Open the  datasets with default value.
    if (argc >= 5)
      datasetId1 = H5Dopen(fileId, v1, H5P_DEFAULT);
    if (argc >= 6)
      datasetId2 = H5Dopen(fileId, v2, H5P_DEFAULT);
    if (argc >= 7)
      datasetId3 = H5Dopen(fileId, v3, H5P_DEFAULT);
    if (argc == 8)  
      datasetId4 = H5Dopen(fileId, v4, H5P_DEFAULT);

    // Create a property list for collective data transfer.
    plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);
    
   
    // Each MPI rank writes its own dataset
    if (rank == 0) 
      status = H5Dread(datasetId1, H5T_NATIVE_INT, mspace, fspace, plistId, data);
    else if (rank == 1)  
      status = H5Dread(datasetId2, H5T_NATIVE_INT, mspace, fspace, plistId, data);
    else if (rank == 2)    
      status = H5Dread(datasetId3, H5T_NATIVE_INT, mspace, fspace, plistId, data);
    else if (rank == 3)    
      status = H5Dread(datasetId4, H5T_NATIVE_INT, mspace, fspace, plistId, data);
    

    if (rank == (size-1))
      printf("Just checking the output %d\n",data[0]); 

    // De-allocate data
    free(data);

    // Close all the HDF5 resources

    H5Dclose(datasetId1);
    if (argc >= 6)
      H5Dclose(datasetId2);
    if (argc >= 7)
      H5Dclose(datasetId3);
    if (argc >= 8)  
      H5Dclose(datasetId4);
    H5Sclose(fspace);
    H5Sclose(mspace);
    H5Fclose(fileId);
    H5Pclose(plistId);
 
    MPI_Finalize();

    return 0;
}     
译