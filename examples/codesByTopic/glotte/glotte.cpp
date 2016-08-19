#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <time.h>
#include "complexDynamics/mrtDynamics.h"
#include "complexDynamics/mrtDynamics.hh"

using namespace plb;
using namespace std;
using namespace plb::descriptors;

//#define DESCRIPTOR descriptors::D3Q19Descriptor
#define DESCRIPTOR descriptors::MRTD3Q19Descriptor

typedef double T;
typedef Array<T,3> Velocity;
typedef MRTdynamics<T,DESCRIPTOR> BackgroundDynamics;
typedef AnechoicMRTdynamics<T,DESCRIPTOR> AnechoicBackgroundDynamics;

#include "acoustics/acoustics3D.h"
using namespace plb_acoustics_3D;


void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter,T rho_in, T rho_out)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");
    
    Box3D Imgbox(0,nx-1, 0, ny-1, util::roundToInt(nz/2-1), util::roundToInt(nz/2-1));

    const T rho0 = 1;
    const T drho = rho0/10;
     
   /* imageWriter.writeScaledGif( createFileName("rhoScaled", iter, 6),
    *computeDensity(lattice, Imgbox),    imSize, imSize );
   */
    imageWriter.writeGif( createFileName("rho", iter, 6),
    *computeDensity(lattice, Imgbox), 
    (T) rho0 - drho/1000000, (T) rho0 + drho/1000000, imSize, imSize);
    // imageWriter.writeGif( createFileName("rho", iter, 6),
    // *computeDensity(lattice, Imgbox),
    // LowBound , UpBound ,imSize, imSize );

    imageWriter.writeScaledGif( createFileName("ux", iter, 6),
    *computeVelocityNorm(lattice, Imgbox), imSize, imSize );

    
    
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
        VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
        vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
        vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);
}
    

int main(int argc, char **argv)
{
    plbInit(&argc, &argv);
    std::string fNameOut = "./tmp/";
    global::directories().setOutputDir(fNameOut);

    //std::string RFName = fNameOut + "density.txt";
    plb_ofstream ResultsFile("./tmp/density.txt");
    plb_ofstream AllSimulationInfo("./tmp/AllSimulationInfo.txt");

////////////////////////////////////VARIAVEIS//////////////////////////////////////

    //// Definição do fluido

    const plint Cs = 330;  
    const T nu = 1.5e-5;
    const T rho0_ar = 1.2;

    ///  Pressão trans-glote

    const T Ptg_cmH20 = 3;
    const T const_tg = 98.0665;
    const T Ptg_Pa = Ptg_cmH20*const_tg;  //[Pa]


    /// Definições do Lattice 

    const plint Res = 12;
    const T lx = 0.04/100;  // [m]
    const T c_LB = 1/ sqrt(3);

    const T dx = lx / Res;
    const T dt = dx * (c_LB/Cs);

    const plint BoudaryLines = 2; // ExtraLinesForBoudaries

    const plint ExtensaoDominio = util::roundUp(((0.3/100) / dx), 2);  // 0.3cm é o comprimento entre entrada e saida laringe

    const plint nx = 432+BoudaryLines+ 3*ExtensaoDominio;
    const plint ny = 258*2+1;  // Rever 
    const plint nz = 8;    
    

    //trocando unidades 
    const T rho0_LB = 1;

    const T rhoTg_LB = (Ptg_Pa / rho0_ar)*((dt*dt)/(dx*dx)) / (c_LB*c_LB)*20 + rho0_LB  ; /// OLHA ESSE 20

    const T nu_LB = nu * dt/(dx*dx);
    const T tau = 0.5 + c_LB*nu_LB;


    // Definições da simulação

    const T omega = 1.9;// tau;

    const T rho_in = rhoTg_LB ;

    const T rho_out = rho0_LB;

    Array<T,3> u0(0, 0, 0);     // Velocidade 


    // Posicionamento da geometria

    T L_in=0.2/100;
    //T L_out=0.1/100;

    T Pos_Bx = 0.7592/100;
    T Pos_By = 0.86/100;

    plint PosX= util::roundToInt(  (L_in+Pos_Bx)/dx ) + BoudaryLines/2;
    plint PosY= util::roundToInt(  (Pos_By)/dx ) ;   

    
    //plint cx = util::roundToInt(nx/2); 
    plint cy = util::roundToInt(ny/2);
    plint cz = util::roundToInt(nz/2);
    //Array<T,3> centerLB(cx , cy, cz);
    Array<T,3> Pos_LB(PosX,PosY-1, +1);

    const plint maxT = 20000;

    // Variaveis para pos-processamento

    plint X1 = 0;  plint X2 = nx-1;
    plint Y1 = cy;  plint Y2 = cy;
    plint Z1 = cz; plint Z2= cz;

    clock_t t;


    // Escrevendo as variaveis da simulação

    AllSimulationInfo <<endl 
    << "Dados da simulação" << endl
    << "Lattice:" << endl << endl 
    << "nx: " << nx << " ny: " << ny << " nz: " << nz << endl
    << "dx: " << dx << " dt: " << dt << endl
    << "tau: " << tau << " omega: " << omega << endl << endl

    << "Tempos: " << endl
    << "Time-steps: " << maxT << endl
    << "Tempo simulado: " << maxT*dt << " s " << endl << endl

    << "Dados do fluido [SI]" << endl
    << "rho: " << rho0_ar << endl
    << "viscosidade cin.: " << nu << endl 
    << "Velocidade do som:" << Cs << endl;

    



/////////////////////////////////CRIAÇAO DO LATTICE////////////////////////////////

    // CRIAÇAO DO LATTICE

    pcout << "Creation of the lattice." << endl;
    pcout << "dx: " << dx << std::endl;
    pcout << "dt: " << dt << std::endl;
    pcout << "tau: " << tau << std::endl;
    pcout << "omega: " << omega << std::endl;

    pcout << std::endl << "Transglottal pressure: " << Ptg_Pa << " Pa" << std::endl;
    pcout << "Transglottal density: " << rho_in << " Lattice units" << std::endl;

    pcout << std::endl << "Number of time-steps: " << maxT << std::endl;
    pcout  << "Total simulated time: " << maxT*dt << " s " << std::endl << std::endl;


 //   MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx,ny,nz, new BGKdynamics<T,DESCRIPTOR>(omega));   /// BGK (nx, ny , nz)
    MultiBlockLattice3D<T, DESCRIPTOR> lattice(nx, ny, nz,  new AnechoicBackgroundDynamics(omega));
    defineDynamics(lattice, lattice.getBoundingBox(), new BackgroundDynamics(omega));
  //  defineDynamics(lattice, lattice.getBoundingBox(), new BackgroundDynamics(omega));

//    SparseBlockStructure3D sparseBlock(nx, ny , nz);
//    sparseBlock.addBlock(Box3D(0, nx-1, 0, ny-1, 1,1),   sparseBlock.nextIncrementalId());
//    sparseBlock.addBlock(Box3D(0, nx-1, 0, ny-1, 2,2),   sparseBlock.nextIncrementalId());
//    sparseBlock.addBlock(Box3D(0, nx-1, 0, ny-1, 3,3),   sparseBlock.nextIncrementalId());
//    sparseBlock.addBlock(Box3D(0, nx-1, 0, ny-1, 3,3),   sparseBlock.nextIncrementalId());

  //  plint envelopeWidth1= 1;

    //MultiBlockLattice3D<T, DESCRIPTOR> lattice ( MultiBlockManagement3D(  sparseBlock ,defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth1 ),
      //  defaultMultiBlockPolicy3D().getBlockCommunicator(),
        //defaultMultiBlockPolicy3D().getCombinedStatistics(),
        //defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
        //new BGKdynamics<T,DESCRIPTOR>(omega)    );

T size_anechoic_buffer = 30; 
plint orientation = 1;
Array<plint, 3 > position_anechoic_wall(nx-1,0,0) ;
plint length_anechoic_wall = ny-1;
plint width_anechoic_wall= nz-1;
T rhoBar_target = 0;
Array<T,3> j_target(0,0,0);

AnechoicBackgroundDynamics *anechoicDynamics = 
new AnechoicBackgroundDynamics(omega);
T delta = 0;
T delta_efective = 30 - delta;
anechoicDynamics->setDelta((T) delta_efective);
anechoicDynamics->setRhoBar_target(rhoBar_target);
//j_target[0] = -j_target[0];  
anechoicDynamics->setJ_target(j_target);
defineDynamics(lattice, Box3D(nx-31, nx-1, 0, ny - 1, 0, nz - 1), anechoicDynamics);


// defineAnechoicMRTWall(nx, ny, nz, lattice,
// size_anechoic_buffer, orientation, omega, 
// position_anechoic_wall, length_anechoic_wall, width_anechoic_wall,
//    rhoBar_target, j_target);


const T rho0 = 1;
const T drho = rho0/10;
initializeAtEquilibrium( lattice, lattice.getBoundingBox(), rho0_LB , u0 );
initializeAtEquilibrium( lattice, Box3D(nx-80,nx-70,cy-10,cy+10,0, nz-1), rho0+drho, u0 );

    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);

    
	// DEFINICACAO DAS C.C.

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    	// ENTRADA E SAIDA

    Box3D  inlet(   0,    0, 2, ny-3, 1,nz-2); // rever
    Box3D outlet(nx-1, nx-1, 2, ny-3, 1,nz-2);

    	// LATERAIS

    Box3D   top(0, nx-1, ny-3, ny-1,    0, nz-1);
    Box3D   bot(0, nx-1,    0,    2,    0, nz-1);
    Box3D front(0, nx-1,    1, ny-2,    0,    0);
    Box3D  back(0, nx-1,    1, ny-2,  nz-1, nz-1);
	
    pcout << "Definition of inlet/outlet." << endl; 

    // boundaryCondition->setPressureConditionOnBlockBoundaries(lattice, inlet, boundary::dirichlet);
    //     setBoundaryDensity(lattice, inlet, rho_in);

    boundaryCondition->setPressureConditionOnBlockBoundaries(lattice, outlet, boundary::dirichlet);
        setBoundaryDensity(lattice, outlet, rho_out);

        // TESTE
    //T Mach = 0.15;
    //T U_c = Mach * c_LB ;
    //T Ve = Res/(ny-2)*U_c;
    //Array<T,3> Vtestes(U_c,0,0);


     //boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, inlet, boundary::dirichlet);
     //boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, outlet, boundary::dirichlet);
     //setBoundaryVelocity(lattice,  inlet, Vtestes );
     //setBoundaryVelocity(lattice, outlet, Vtestes );

    pcout << "Generating outter boundary conditions." << std::endl;

    //defineDynamics(lattice, top,   new BounceBack<T,DESCRIPTOR>());
    //defineDynamics(lattice, bot,   new BounceBack<T,DESCRIPTOR>());
    //defineDynamics(lattice, front, new BounceBack<T,DESCRIPTOR>());
    //defineDynamics(lattice,  back, new BounceBack<T,DESCRIPTOR>());


    defineDynamics(lattice, top,   new NoDynamics<T,DESCRIPTOR>());
    defineDynamics(lattice, bot,   new NoDynamics<T,DESCRIPTOR>());
    
    defineDynamics(lattice, front, new NoDynamics<T,DESCRIPTOR>());
    defineDynamics(lattice, back,  new NoDynamics<T,DESCRIPTOR>());

    //boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, top  , boundary::freeslip);
    //boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, bot  , boundary::freeslip);    
    //boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, front, boundary::freeslip);
    //boundaryCondition->setVelocityConditionOnBlockBoundaries(lattice, back , boundary::freeslip);
    
        //setBoundaryVelocity(lattice,   top, u0);
        //setBoundaryVelocity(lattice,   bot, u0);
        //setBoundaryVelocity(lattice, front, u0);
        //setBoundaryVelocity(lattice,  back, u0);

        // Condições iniciais
        
	pcout << "Initialization of rho and u." << endl;
    

    lattice.initialize();
    //delete boundaryCondition;

////////////////////////////////////GEOMETRIAS/////////////////////////////////////


    // Ler a geometria
   
    plint xDirection = 0;
    plint borderWidth = 1;      
                                
    plint margin = 1;           
    plint blockSize = 0;   

    plint extendedEnvelopeWidth = 2;   // Extrapolated off-lattice BCs.
    const int flowType = voxelFlag::outside;

    bool useAllDirections=true;
    BoundaryProfiles3D<T,Velocity> profiles;
    profiles.setWallProfile(new NoSlipProfile3D<T>);



    pcout << std::endl << "Reading STL data for the obstacle geometry 1." << std::endl;   
    TriangleSet<T> triangleSet("glote_up_1_50mm.STL", DBL);

    Array<T,3> obstacleCenter(0,0,0);
    triangleSet.scale(1/dx /100);  // Essa divisao por 100 é devido a diferença de unidade no Solid 
    triangleSet.translate(Pos_LB);

    triangleSet.writeBinarySTL(fNameOut+"/Prega_Up.stl");

    //triangleSet.refineRecursively(dx, 3);
   
    DEFscaledMesh<T> defMesh(triangleSet, 0, xDirection, margin, Dot3D(0, 0, 0));
    TriangleBoundary3D<T> boundary(defMesh);

        // Voxelizando o dominio

    pcout << "Voxelizing the domain 1" << std::endl;
  
    VoxelizedDomain3D<T> voxelizedDomain ( boundary, flowType, lattice.getBoundingBox(), borderWidth, extendedEnvelopeWidth, blockSize );         
    defineDynamics(lattice, voxelizedDomain.getVoxelMatrix(), lattice.getBoundingBox(), new NoDynamics<T,DESCRIPTOR>() , voxelFlag::inside);  // new BounceBack<T,DESCRIPTOR>()

    // Condição de contorno não alinhada

    pcout << "Generating non-alligned boundary conditions 1." << std::endl;

    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity> *OffboundaryCondition;
    OffLatticeModel3D<T,Velocity>* offLatticeModel=0;
        
    offLatticeModel = new GuoOffLatticeModel3D<T,DESCRIPTOR>( new TriangleFlowShape3D<T,Array<T,3> >(voxelizedDomain.getBoundary(), profiles), flowType, useAllDirections );
    OffboundaryCondition = new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>( offLatticeModel, voxelizedDomain, lattice);

    OffboundaryCondition->insert();
    
    pcout << "Done." << std::endl << std::endl;


    pcout << std::endl << "Reading STL data for the obstacle geometry 2." << std::endl;   
    TriangleSet<T> triangleSet2("glote_down_1_50mm.STL", DBL);

    triangleSet2.scale(1/dx /100);  // Essa divisao por 100 é devido a diferença de unidade no Solid 
    triangleSet2.translate(Pos_LB);

    triangleSet2.writeBinarySTL(fNameOut+"/Prega_Down.stl");

    //triangleSet.refineRecursively(dx, 3);
   
    DEFscaledMesh<T> defMesh2(triangleSet2, 0, xDirection, margin, Dot3D(0, 0, 0));
    TriangleBoundary3D<T> boundary2(defMesh2);

        // Voxelizando o dominio

    pcout << "Voxelizing the domain 2" << std::endl;
  
    VoxelizedDomain3D<T> voxelizedDomain2 ( boundary2, flowType, lattice.getBoundingBox(), borderWidth, extendedEnvelopeWidth, blockSize );         
    defineDynamics(lattice, voxelizedDomain2.getVoxelMatrix(), lattice.getBoundingBox(), new NoDynamics<T,DESCRIPTOR>() , voxelFlag::inside);  // new NoDynamics<T,DESCRIPTOR>()

    // Condição de contorno não alinhada

    pcout << "Generating non-alligned boundary conditions 2." << std::endl;

    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity> *OffboundaryCondition2;
    OffLatticeModel3D<T,Velocity>* offLatticeModel2=0;
        
    offLatticeModel2 = new GuoOffLatticeModel3D<T,DESCRIPTOR>( new TriangleFlowShape3D<T,Array<T,3> >(voxelizedDomain2.getBoundary(), profiles), flowType, useAllDirections );
    OffboundaryCondition2 = new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>( offLatticeModel2, voxelizedDomain2, lattice);

    OffboundaryCondition2->insert();
    
    pcout << "Done." << std::endl << std::endl;

    pcout << "Parallel info:" << std::endl ;
    MultiBlockManagement3D sparseBlockManagement(lattice.getMultiBlockManagement());
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    lattice.toggleInternalStatistics(false);


        
//////////////////////////////////////MAIN LOOP////////////////////////////////////

    // INICIO DA SIMULAÇAO PER-SE

    pcout << "Simulation begins" << endl;
    plint iT=0;
    t = clock();
    
    for (iT=0;iT<maxT; ++iT)
    {
                    
        if (iT % 20 == 0) {
            pcout << "Iteration " << iT << endl;
        }

        if (iT % 20 == 0 ) {
            pcout << "Writing Gif"<< endl;
            writeGifs(lattice,iT,rho_in,rho_out);
        }

        if(iT % 250 == 0){
            pcout << "Writing VTK"<< endl;
        	//writeVTK(lattice, iT);
        }

        lattice.collideAndStream();

    }

    t = (clock() - t)/CLOCKS_PER_SEC;
    AllSimulationInfo << endl << "Main-loop Execution time: " << t << endl << endl;

    pcout << "End of simulation at iteration " << iT << endl;


//////////////////////////////////////POSPROCESSAMENTO/////////////////////////////

    // POS-PROCESSAMENTO

    pcout << endl << " Creating files for post-processing" << endl ;

    pcout << " Writing results along: ( " << X1 << ", "<< X2 << ", "<< Y1 << ", "<< Y2 << ", "<< Z1 << ", "<< Z2 << ") " << endl ;
    Box3D ResLine( X1 , X2, Y1, Y2 , Z1 , Z2 );
    ResultsFile <<  setprecision(10) << *computeDensity(lattice,ResLine) << endl;
    ResultsFile.close();
    AllSimulationInfo.close(); 




    pcout << "End of Program" << endl;

}

