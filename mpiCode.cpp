#include <cstdlib>
#include <mpi.h>
#include <vector>
#include <unistd.h>
#include <string>
#include <fstream>
#include <math.h>
#include <matrix.h>

using std::vector;

// function Prototypes 
void FillFieldGaussian( Matrix<double> &F, vector<double> xCoord,vector<double> yCoord, int Nghost);
void ExchangeIBorders( Matrix<double> &F,int coordpXaxis, int ni, int nj, int NpI, int Nghost, int lowIneighbor, int highIneighbor);
void ExchangeJBorders( Matrix<double> &F,int coordpYaxis, int ni, int nj, int NpJ, int Nghost, int lowJneighbor, int highJneighbor);
void PrintField( Matrix<double> F,int rank);

// MPI Send and Recieve Tags 
enum {F_GHOST_TAG, I_GHOST_TAG_0, I_GHOST_TAG_1, J_GHOST_TAG_0, J_GHOST_TAG_1};

int static TRUE = 1; int static FALSE = 0; //#todo 
int static XAXIS = 0; int static YAXIS = 1; int static ZAXIS = 2; //enum

//global MPI variables (I KNOW!!!)
MPI_Status status;          // mpi: store status of a MPI_Recv
MPI_Request request;        // mpi: capture request of a MPI_Isend
MPI_Datatype colType;  
int main(int argc, char *argv[]) {
  int rank;                   // mpi: process id number
  int nProcesses;             // mpi: number of total processess 
  MPI_Init(&argc, &argv);		      /* initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD, &nProcesses); /* store the number of processes */  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);	      /* store the rank */

  // The information each process needs (local info)
  vector<double> xCoord; // holds local xCoordinates
  vector<double> yCoord; // holds local xCoordinates  

  int Nghost;            // defines the number of ghost cells (ONLY WORKS FOR ONE!!1 #todo)
  int ni;                // holds the number of i-points calculated by local process, 
  double dx;
  int nj;                // holds the number of i-points calculated by local process, 
  double dy;
  int NpI, NpJ;

  // READ IN from parameters.input #todo
  int static NDIMS = 2; 		// number of dimensions
  int static Ni = 8; 		// number of global grid points in i-direction
  int static Nj = 6; 		// number of global grid points in j-direction  
  Nghost = 1;	        // number of points in ghost cells  #todo only works for 1!!!
  double xmin = 0;	// domain min
  double xmax = 7;	// domain max
  double ymin = 0;	// domain min
  double ymax = 24;	// domain max
  
  // set the number of processes in each dimension, recieved as input arguements
  if (NDIMS >0){
    NpI = atoi(argv[1]);
    if (rank == 0) std::cout << "NpI = " << NpI << "\n";
    if (NDIMS >1){
      NpJ = atoi(argv[2]);
      if (rank == 0) std::cout << "NpJ = " << NpJ << "\n";
    } 
  }
  // Check some errors
  if (NpI*NpJ != nProcesses){
    if (rank == 0) std::cout << "ERROR: Nx * Ny must equal Number of Processes \n";
    return -1;
  }
  if (Ni % NpI != 0){
    if (rank == 0) std::cout << "ERROR: Nx divisible by NpI \n";
    return -1;
  }
  if (Ni % NpJ != 0){
    if (rank == 0) std::cout << "ERROR: Ny divisible by NpJ \n";
    return -1;
  }

  ni = (Ni / (NpI)); //expected size of the parallel chunk  
  //Error check: size of the calculations should be larger than the ghost-cell region
  if (Nghost >= ni) {
    if (rank == 0) std::cout << "ERROR: Please set the number of Nghost << Granularity. \n  " ;
    return -1;
  }    
  nj = (Nj / (NpJ)); //expected size of the parallel chunk  
  //Error check: size of the calculations should be larger than the ghost-cell region
  if (Nghost >= nj) {
    if (rank == 0) std::cout << "ERROR: Please set the number of Nghost << Granularity. \n  " ;
    return -1;
  }    

  //caculated from input
  dx = (xmax -xmin)/(Ni-1);
  //Global Info: created and deleted
  // created by all, only written by rank ==0 
  std::vector<double> globalXcoord(Ni);
  for(int i=0; i < Ni; i++){
    globalXcoord[i] = xmin + i*dx;
  }
  //caculated from input
  dy = (ymax -ymin)/(Nj-1);
  //Global Info: created and deleted
  // created by all, only written by rank ==0 
  std::vector<double> globalYcoord(Nj);
  for(int i=0; i < Nj; i++){
    globalYcoord[i] = ymin + i*dy;
  }
  double dt = 1
  // comm cart
  int neighbor[NDIMS][2],     /* ranks of neighbors in each coordinate */
      np[NDIMS],              /* number of processors in each coordinate */
      coordp[NDIMS];          /* coordinate of the processor in the cpu grid */

   MPI_Type_vector(ni,		/* # column elements */
		   1,			/* 1 column only */
		   (nj+2*Nghost),	/* skip 20 elements */
		   MPI_DOUBLE,		/* elements are float */
		   &colType);		/* MPI derived datatype */

  // MPI_Type_vector(6,		/* # column elements */
  // 		  1,			/* 1 column only */
  // 		  14/6,	/* skip 20 elements */
  // 		  MPI_DOUBLE,		/* elements are float */
  // 		  &colType);		/* MPI derived datatype */

   MPI_Type_commit(&colType);

  MPI_Comm MPI_COMM_CART;
  int periodic[NDIMS], cartps[NDIMS], reorder;
  reorder=FALSE;  //reorder the processes? 
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  np[XAXIS]=NpI; np[YAXIS] = NpJ;
  periodic[XAXIS]=TRUE; periodic[YAXIS]=TRUE; //periodic?
  MPI_Cart_create(MPI_COMM_WORLD, NDIMS, np, periodic, reorder, &MPI_COMM_CART);
  /* get own position */
  MPI_Cart_get(MPI_COMM_CART, NDIMS, np, periodic, coordp);
  /* get own rank */
  MPI_Cart_rank(MPI_COMM_CART, coordp, &rank);
  
  //set the low and high neighbors in each dimension. If there is no neighbor in a certian direction, the value will be set to -2.
  MPI_Cart_shift(MPI_COMM_CART, XAXIS, 1, &neighbor[XAXIS][0], &neighbor[XAXIS][1]);
  MPI_Cart_shift(MPI_COMM_CART, YAXIS, 1, &neighbor[YAXIS][0], &neighbor[YAXIS][1]);
  // sleep(1*rank);
  // std::cout  << "rank I: " << rank<< " lowIneighbor: " << neighbor[XAXIS][0] << " highIneighbor: " << neighbor[XAXIS][1] <<" \n";
  // std::cout  << "rank J: " << rank<< " lowJneighbor: " << neighbor[YAXIS][0] << " highJneighbor: " << neighbor[YAXIS][1] <<" \n";
  // std::cout  << "rank: " << rank<< " pocX: " << coordp[XAXIS] << " porcY: " << coordp[YAXIS] <<" \n";    
  
  
  ///////////* Master initializes work (work that must only be done once)*///////////
  if (rank == 0) {  
    std::cout << "Starting an MPI parallel Lapalcian.  \n  " ;
    std::cout << "Reading parameters.input  (hardcoded 1D) #todo.  \n  " ;
    std::cout << "Setting up the global uniform grid.  \n  " ;    
    std::cout << "Write out global xCoord to file. \n" ;

    std::ofstream fsXcoord( "xcoordinates.out" );
    for(vector<double>::const_iterator i = globalXcoord.begin(); i != globalXcoord.end(); ++i) {
      fsXcoord << *i << '\n';
    }
    std::ofstream fsYcoord( "ycoordinates.out" );
    for(vector<double>::const_iterator i = globalYcoord.begin(); i != globalYcoord.end(); ++i) {
      fsYcoord << *i << '\n';
    }
  }
  
  int iStart = coordp[XAXIS] * ni;      //set the start index (in global grid)
  int iEnd = iStart + ni;
  int jStart = coordp[YAXIS] * nj;      //set the start index (in global grid)
  int jEnd = jStart + nj;
        
  //fill local x-coord array
  xCoord.reserve(ni); // holds local xCoordinates
  for( int i = iStart; i<iEnd; i++){
    xCoord.push_back(globalXcoord[i]);
  }

  //fill local x-coord array
  yCoord.reserve(nj); // holds local xCoordinates
  for( int i = jStart; i<jEnd; i++){
    yCoord.push_back(globalYcoord[i]);
  }
  
  //#todo delete globalXcoord at this point
  
  // at this point every process has its xCoord, ni, and Nghost  (mpi setup complete)
#ifdef DEBUG
  std::cout << "rank: " << rank << " ni: " << ni << " Nghost: " << Nghost << " \n" ;
  std::cout << "Print x-coords. \n" ;
  for(int i=0;  i < ni; i++){
    std::cout << rank<< ": "<< xCoord[i] << "\n" ;          
    }
  std::cout << "rank: " << rank << " nj: " << nj << " Nghost: " << Nghost << " \n" ;
  std::cout << "Print y-coords. \n" ;
  for(int i=0;  i < nj; i++){
    std::cout << rank<< ": "<< yCoord[i] << "\n" ;          
    }
  
  // sleep(2*rank); //seperate the output of the processors in the terminal
#endif
  std::cout << "initialize field \n";
  Matrix<double> field = Matrix<double>(ni+2*Nghost,nj+2*Nghost);
  Matrix<double> fieldNew = Matrix<double>(ni+2*Nghost,nj+2*Nghost); 
  
  ////////////////////fill field ////////////////
  FillFieldGaussian(field, xCoord, yCoord, Nghost);

  std::cout << "save intitial field \n";
  //save data
  std::string initialFileName = "inField"+std::to_string(rank)+".out";   
  std::ofstream initialF(initialFileName);

  for (int j = 0; j < field.cols(); j++)  {
    initialF << '\n';
    for (int i = 0; i < field.rows(); i++) 
      initialF << field(i,j) << ' ';
    }

  for(int T = 0; T < 1000; T++) {
    /////////send boundaries (ghost cells)////////
    ExchangeIBorders(field, coordp[XAXIS], ni, nj,  NpI,  Nghost, neighbor[XAXIS][0], neighbor[XAXIS][1]);
    ExchangeJBorders(field, coordp[YAXIS], ni, nj,  NpJ,  Nghost, neighbor[YAXIS][0], neighbor[YAXIS][1]); 

    //   ////////////////Laplacian//////////////////////
    for (int i = Nghost; i < ni+Nghost; i++){
      for (int j = Nghost; j < nj+Nghost; j++){    
	fieldNew(i,j) =field(i,j)+ dt*((field(i-1,j) -2*field(i,j) + field(i+1,j))/(2*dx)
				     + (field(i,j-1) - 2*field(i,j) + field(i,j+1))/(2*dy));
      }
    }

  }
  std::cout << "save Final field \n";
  ///////////////save data to file/////////////
  std::string finalFileName = "finalField"+std::to_string(rank)+".out";   
  std::ofstream finalF( finalFileName);
  for (int j = 0; j < fieldNew.cols(); j++){ 
    finalF << '\n';
    for (int i = 0; i < fieldNew.rows(); i++) 
      finalF << fieldNew(i,j) << ' ';
  }  

  // field = fieldNew;
  MPI_Finalize(); //finalize MPI operations
  return 0;
}

void ExchangeHighJ(Matrix<double> &F, int ni, int nj, int Nghost, int lowJneighbor, int highJneighbor){
    MPI_Isend(&F(Nghost,    nj   ), 1, colType, highJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &request);
    MPI_Recv( &F(Nghost,    nj+Nghost    ), 1, colType, highJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &status);      

  // MPI_Sendrecv(&F(Nghost,   nj    ), 1, colType, highJneighbor, J_GHOST_TAG,
  // 	       &F(Nghost,nj+Nghost), 1, colType, highJneighbor, J_GHOST_TAG, MPI_COMM_WORLD, &status);
}
void ExchangeLowJ( Matrix<double> &F,int ni, int nj, int Nghost, int lowJneighbor, int highJneighbor){
    MPI_Recv( &F(Nghost,    0    ), 1, colType, highJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &status);
    MPI_Isend(&F(Nghost,    Nghost   ), 1, colType, highJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &request);    

 //  MPI_Sendrecv(&F(Nghost, Nghost), 1, colType, lowJneighbor, J_GHOST_TAG,
 // 	       &F(Nghost,   0   ), 1, colType, lowJneighbor, J_GHOST_TAG, MPI_COMM_WORLD, &status);


}


void ExchangeJBorders( Matrix<double> &F, int coordpYaxis, int ni, int nj, int NpJ, int Nghost, int lowJneighbor, int highJneighbor){
  if (NpJ ==1){
    MPI_Isend(&F(Nghost,    nj   ), 1, colType, highJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &request);
    MPI_Recv( &F(Nghost,    0    ), 1, colType, lowJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &status);
    MPI_Isend(&F(Nghost,    Nghost   ), 1, colType, highJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &request);
    MPI_Recv( &F(Nghost,    nj+Nghost    ), 1, colType, lowJneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &status);  
    
  }
 else if (coordpYaxis % 2 == 0){
    ExchangeLowJ( F, ni, nj, Nghost, lowJneighbor, highJneighbor);       
    ExchangeHighJ(F, ni, nj, Nghost, lowJneighbor, highJneighbor);      
  }
  else{
    ExchangeHighJ(F, ni, nj, Nghost, lowJneighbor, highJneighbor);      
    ExchangeLowJ( F, ni, nj, Nghost, lowJneighbor, highJneighbor);           
    }
}


void ExchangeHighI(Matrix<double> &F, int ni, int nj, int Nghost, int lowIneighbor, int highIneighbor){
  // MPI_Isend(&F(ni   ,    Nghost), Nghost*nj, MPI_DOUBLE, highIneighbor ,F_GHOST_TAG, MPI_COMM_WORLD,&request);
  // MPI_Recv(&F(ni+Nghost, Nghost), Nghost*nj, MPI_DOUBLE, highIneighbor ,F_GHOST_TAG,
  MPI_COMM_WORLD, &status);          
  //above done in one command
  MPI_Sendrecv(&F(  ni   ,Nghost), Nghost*nj, MPI_DOUBLE, highIneighbor , 2,
  	       &F(ni+Nghost,Nghost), Nghost*nj, MPI_DOUBLE, highIneighbor,2, MPI_COMM_WORLD, &status);          

}
void ExchangeLowI( Matrix<double> &F, int ni, int nj, int Nghost, int lowIneighbor, int highIneighbor){
  // MPI_Recv( &F(   0  , Nghost), Nghost*nj, MPI_DOUBLE, lowIneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &status);
  // MPI_Isend(&F(Nghost, Nghost), Nghost*nj, MPI_DOUBLE, lowIneighbor, F_GHOST_TAG, MPI_COMM_WORLD,&request);
  //above done in one command  
  MPI_Sendrecv(&F(Nghost, Nghost)  , Nghost*nj, MPI_DOUBLE, lowIneighbor, 2,
  	       &F(  0    ,  Nghost), Nghost*nj, MPI_DOUBLE, lowIneighbor, 2, MPI_COMM_WORLD, &status);	            		       

  }


void ExchangeIBorders( Matrix<double> &F, int coordpXaxis, int ni, int nj, int NpJ, int Nghost, int lowIneighbor, int highIneighbor){
  if (NpJ == 1){
    MPI_Isend(&F(ni   ,    Nghost), Nghost*nj, MPI_DOUBLE, highIneighbor ,F_GHOST_TAG, MPI_COMM_WORLD,&request);
    MPI_Recv( &F(   0  , Nghost), Nghost*nj, MPI_DOUBLE, lowIneighbor, F_GHOST_TAG, MPI_COMM_WORLD, &status);    
    MPI_Isend(&F(Nghost, Nghost), Nghost*nj, MPI_DOUBLE, lowIneighbor, F_GHOST_TAG, MPI_COMM_WORLD,&request);  
    MPI_Recv(&F(ni+Nghost, Nghost), Nghost*nj, MPI_DOUBLE, highIneighbor ,F_GHOST_TAG, MPI_COMM_WORLD, &status);
  }
  else if (coordpXaxis % 2 == 0){
    ExchangeLowI( F, ni, nj, Nghost,  lowIneighbor, highIneighbor);      
    ExchangeHighI(F, ni, nj, Nghost,  lowIneighbor, highIneighbor);      
  }
  else{
    ExchangeHighI(F, ni, nj, Nghost,  lowIneighbor, highIneighbor);      
    ExchangeLowI( F, ni, nj, Nghost,  lowIneighbor, highIneighbor);      
    }
}

void FillFieldGaussian(Matrix<double> &F,vector<double> xCoord, vector<double> yCoord, int Nghost){
  /* fix wonky indexing; possibly modify the Matrix class to hold Nghost?  #todo */  
  for(int i = Nghost; i < (F.rows()-Nghost);i++){
    for(int j = Nghost; j < (F.cols()-Nghost);j++){
      // F(i,j) = exp(- sqrt(xCoord[i-Nghost]*xCoord[i-Nghost] +yCoord[j-Nghost]*yCoord[j-Nghost]) /10.0) ;
      F(i,j) = xCoord[i-Nghost]+yCoord[j-Nghost];
      // F(i,j) = (i-Nghost)+(j-Nghost)*8;      
      }
  }
}
