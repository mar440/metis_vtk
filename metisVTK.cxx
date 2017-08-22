//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                                                     +
//     DECOMPOSITION OF 3D MESH USING METIS                                            +
//     ===================================                                             +
//     author:     A. Markopuolos                                                      +
//     revison:    20170822                                                            +
//                                                                                     +
//     LAUNCHING:                                                                      +
//                                                                                     +
//         $ ./metisVTK path_to_vtu_file X                                             +
//                                                                                     +
//     with parameters + path_to_vtu_file - path to the vtu file with mesh,            +
//                     + X - required number of subdomains.                            +
//                                                                                     +
//     NOTE:                                                                           +
//     If mesh contains disjoned parts, it may happen,                                 +
//     required number of subdomains will be not satisfied (due to load balancig).     +
//                                                                                     +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vector>
#include <math.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include <string>


// To extract a sub part
#include <vtkThreshold.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>
#include <vtkInformation.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkAppendFilter.h>
#include <set>
#include <stdexcept>
#include "metis.h"
#include <vtkConnectivityFilter.h>
//#include <vtkIdTypeArray.h>
#include <iostream>

using namespace std;
bool ASCII_OR_BINARY_VTU = true;


void printVTK(vtkSmartPointer<vtkUnstructuredGrid> _unstructuredGrid,string fname){
    cout <<"vtkVersion: " << vtkVersion::GetVTKMajorVersion() << "\n";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fname.c_str());
    if (ASCII_OR_BINARY_VTU){
        writer->SetDataModeToAscii();
    }
    else{
        writer->SetDataModeToBinary();
    }
    writer->SetInputData(_unstructuredGrid);
    writer->Write();
}


void decomposeRegion(vtkSmartPointer<vtkUnstructuredGrid> mesh_local, vtkSmartPointer<vtkDataArray> globalPartitionId, int _nparts, int prevPartitionId){

    int nPointsLocal = mesh_local->GetNumberOfPoints();
    int nCellsLocal = mesh_local->GetNumberOfCells();

    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    int nPi;
    int * eptr = new int[nCellsLocal + 1];
    eptr[0] = 0;
    for (int i = 0 ; i < nCellsLocal; i++){
        mesh_local->GetCell(i,cell);
        nPi = cell->GetNumberOfPoints();
        eptr[i+1] = eptr[i] + nPi;
    }
    int * eind = new int[eptr[nCellsLocal]];
    int cnt = 0;
    int CellDim;
    for (int i = 0 ; i < nCellsLocal; i++){
        mesh_local->GetCell(i,cell);
        nPi = cell->GetNumberOfPoints();
        for (int j = 0; j < nPi; j++){
            eind[cnt] = cell->GetPointId(j);
            cnt ++;
        }
        if (CellDim <cell->GetCellDimension()){
            CellDim  = cell->GetCellDimension();
        }
    }

    int ncommon = 2; /* ncommon = 2 for all other types ... */
    if (CellDim == 2){
        ncommon = 2;
    }
    else if (CellDim == 3){
        ncommon = 3;
    }

    int options[METIS_NOPTIONS];
    cout << " ----------   METIS_NOPTIONS " << METIS_NOPTIONS << endl;
    options[METIS_OPTION_PTYPE    ] = METIS_PTYPE_RB;    // multilevel recursive bisectioning
    options[METIS_OPTION_OBJTYPE  ] = METIS_OBJTYPE_CUT; // edge-cut minimization
    options[METIS_OPTION_CTYPE    ] = METIS_CTYPE_RM;    // random matching
    options[METIS_OPTION_IPTYPE   ] = METIS_IPTYPE_GROW; // grows a bisction using a greedy strategy
    options[METIS_OPTION_RTYPE    ] = METIS_RTYPE_FM;    // FM-based cut refinement
    options[METIS_OPTION_NCUTS    ] = 1 ;//
    options[METIS_OPTION_NITER    ] = 10;/* Default value */
    options[METIS_OPTION_SEED     ] = -1;/* Seed of random algo */
    options[METIS_OPTION_UFACTOR  ] = 1;
    options[METIS_OPTION_NUMBERING] = 0; // C-style numbering
    options[METIS_OPTION_DBGLVL   ] = METIS_DBG_INFO;
    options[METIS_OPTION_CONTIG   ] = 1;

    int nparts = 1;
    int objval;


    if (_nparts > 0 ){
       nparts = _nparts;
    }

    int *epart = new int [nCellsLocal];
    int *npart = new int [eptr[nCellsLocal]];

    METIS_PartMeshDual(&nCellsLocal,    // number of elements in the mesh       Y
                       &nPointsLocal,   //                                      Y
                       eptr,            //                                      Y
                       eind,            //                                      Y
                       (int *)NULL,     // vwgt                                 Y
                       (int *)NULL,     // vsize                                Y
                       &ncommon,        //                                      Y
                       &nparts,         //                                      N
                       (real_t*)NULL,   // tpwgts                               Y
                       options,         //                                      Y
                       &objval,         //                                      Y
                       epart,           //                                      N
                       npart);          //                                      Y


    double tuple[] = {0};
    int _i;

    vtkDataArray * localPartitionId;
    localPartitionId = mesh_local->GetCellData()->GetArray("PartitionId");


    for (int i = 0 ; i < nCellsLocal; i++){
        /* global PartitionId */
        tuple[0] = epart[i] + prevPartitionId;
        /* first: recovering the position of cell in original 'mesh_local' */
        _i = localPartitionId->GetTuple1(i);
        /* second: replacing by localPartitionId + (number of partitions of all previous Regions) */
        globalPartitionId->SetTuple(_i,tuple);
    }

    delete [] epart;
    delete [] npart;
    delete [] eptr;
    delete [] eind;

}


int main(int argc, char *argv[])
{
    if (argc != 3){
        cout << "!!!  Number of arguments is incorrect !!!\n" << endl;
        string str0 = "/path/to/mesh.vtu";
        cout << "To decompose the mesh (stored in "<< str0 <<") into, e.g., " <<
                "3 subdomains, call:\n\n\n \t\t\t./metisVTK "<<str0<< " 3\n\n\n";
        return 0;
    }

    string filename = argv[1];
    string _nparts_str = argv[2];

    int _nparts = stoi(_nparts_str);

    cout << "numb. of inp.             " << argc << endl;
    cout << "file name to be readed is " << filename << endl;
    cout << "nparts from command line: " << _nparts << endl;

    int nparts;
    if (_nparts > 0 ){
       nparts = _nparts;
    }


    //read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> meshGlobal =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    meshGlobal->SetFileName(filename.c_str());
    meshGlobal->Update();

    int _i, nCellArr = meshGlobal->GetNumberOfCellArrays();
    string str_i;
    bool flagPartitionId = false;
    for (_i = 0; _i < nCellArr ; _i++) {
        str_i = meshGlobal->GetCellArrayName(_i);
        if (str_i.compare("PartitionId") == 0){
            flagPartitionId = true;
            break;
        }
    }
    int nPointsGlobal = meshGlobal->GetNumberOfPoints();
    int nCellsGlobal= meshGlobal->GetNumberOfCells();
    cout << "meshGlobal->GetNumberOfPoints() = " << nPointsGlobal << endl;
    cout << "meshGlobal->GetNumberOfCells()" << nCellsGlobal<< endl;




    //vtkDataArray * PartitionId;
    //vtkSmartPointer< vtkDataArray > PartitionId;
    vtkSmartPointer< vtkIntArray > PartitionId;


    if (flagPartitionId){
//        PartitionId = meshGlobal->GetOutput()->GetCellData()->GetArray("PartitionId");
        PartitionId = vtkIntArray::SafeDownCast(meshGlobal->GetOutput()->GetCellData()->GetArray("PartitionId"));
    }
    else{
//        PartitionId = vtkIntArray::New();
        PartitionId = vtkIntArray::New();
        meshGlobal->GetOutput()->GetCellData()->AddArray(PartitionId);
        PartitionId->SetName("PartitionId");
        PartitionId->SetNumberOfComponents(1);
        PartitionId->SetNumberOfTuples(nCellsGlobal);
    }

    /* Partition is firstly used to store original order of cells */
    double tuple[] = {0};
    for (int i = 0 ; i < nCellsGlobal; i++){
        tuple[0] = i;
        PartitionId->SetTuple(i,tuple);
    }



    vtkSmartPointer<vtkConnectivityFilter> connectivityFilter =
    vtkSmartPointer<vtkConnectivityFilter>::New();
    //connectivityFilter->SetInputConnection(meshGlobal->GetOutputPort());
    connectivityFilter->SetInputData(meshGlobal->GetOutput());
    connectivityFilter->SetExtractionModeToAllRegions();
    connectivityFilter->ColorRegionsOn();
    connectivityFilter->Update();
    cout << " MEMORY   after del: " << connectivityFilter->GetOutput()->GetActualMemorySize() << endl;

    int nRegions = connectivityFilter->GetNumberOfExtractedRegions();




    /* First: get number of cells per each Region */
    vector < int > nCellsPerRegion;
    nCellsPerRegion.resize(nRegions);

    vector < double > nPartsPerRegion;
    nPartsPerRegion.resize(nRegions);

    double nAverCellsPerRegion = double(nCellsGlobal) / nparts;
        cout <<"aver. # el. per sub = " << nAverCellsPerRegion  << endl;

        cout <<"=========================================="<< endl;
        cout <<"Global numb.of cells: "<<nCellsGlobal << endl;
        cout <<"=========================================="<< endl;
        cout <<"given numb. of parts = " << nparts << endl;
    double __nparts;

    for (int i = 0; i < nRegions; i ++) {
        vtkSmartPointer <vtkUnstructuredGrid> omega_i = vtkSmartPointer< vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
        threshold->SetInputData(connectivityFilter->GetOutput());
        threshold->ThresholdBetween((double)i, (double)i);
        threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
        threshold->Update();
        omega_i->ShallowCopy(threshold->GetOutput());

        nCellsPerRegion[i] = omega_i->GetNumberOfCells();
        __nparts = double( nCellsPerRegion[i]) / nAverCellsPerRegion;
        nPartsPerRegion[i] = round(__nparts);
        if (nPartsPerRegion[i] == 0)
            nPartsPerRegion[i] = 1;

        cout <<".........................................."     << endl;
        cout <<"nPartsPerRegion  =    " << __nparts             << endl;
        cout <<"nPartsPerRegion  =    " << nPartsPerRegion[i]   << endl;
        cout <<"nCellsPerRegion  =    " << nCellsPerRegion[i]   << endl;
    }

        cout <<"=========================================="<< endl;



    int prevPartitionId = 0;
    for (int i = 0; i < nRegions; i ++) {
        vtkSmartPointer <vtkUnstructuredGrid> omega_i = vtkSmartPointer< vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
        threshold->SetInputData(connectivityFilter->GetOutput());
        threshold->ThresholdBetween((double)i, (double)i);
        threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
        threshold->Update();
        omega_i->ShallowCopy(threshold->GetOutput());
        cout << "omega_i->GetActualMemorySize("<<i<<") = " << omega_i->GetActualMemorySize() << endl;

        decomposeRegion(omega_i, PartitionId, nPartsPerRegion[i], prevPartitionId);
        prevPartitionId += nPartsPerRegion[i];

        //  string fname2 = "dmpFls/Omega_"+to_string(i)+".vtu";
        //  printVTK(omega_i,fname2);
    }

    int real_nparts = prevPartitionId;

    /* Extraction of original file name from the 'path'*/
    stringstream ss1(filename);
    vector <string> result;

    while (ss1.good() )
    {
        string substr;
        getline (ss1, substr, '/');
        result.push_back( substr );
    }

    stringstream ss2(result[result.size() - 1]);
    result.clear();
    result.shrink_to_fit();

    while (ss2.good() )
    {
        string substr;
        getline (ss2, substr, '.');
        result.push_back( substr );
    }
    cout << result.size() << endl;
    cout << result[0]<< endl;
  /* ------------------ */
    string fname = result[0] + "_" + to_string(real_nparts) + "_subs.vtu";

    vtkSmartPointer<vtkUnstructuredGrid> ModelMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ModelMesh->ShallowCopy(meshGlobal->GetOutput());
    printVTK(ModelMesh,fname);


#if 0
    int  nPointArr = meshGlobal->GetNumberOfPointArrays();
    str_i;
    bool flagdisplace = false;
    for (_i = 0; _i < nPointArr ; _i++) {
        str_i = meshGlobal->GetPointArrayName(_i);
        if (str_i.compare("displace") == 0){
            flagdisplace = true;
            break;
        }
    }
    vtkSmartPointer< vtkDataArray > displace;
    if (flagdisplace){
        displace = meshGlobal->GetOutput()->GetPointData()->GetArray("displace");
        ofstream myfile;
        string fname1 = result[0] + "_displace.txt";
        myfile.open (fname1);
        string dataType = displace->GetDataTypeAsString();
        cout << " ++++++++++++++++++++++  " << dataType << endl;
        cout << " ++++++++++++++++++++++  " << displace->GetDataSize()<< endl;
        cout << " ++++++++++++++++++++++  " << displace->GetNumberOfComponents()<< endl;
        cout << " ++++++++++++++++++++++  " << displace->GetSize()<< endl;
        double  u_i[3];
        double  *pu_i;
        pu_i = u_i;
        for (int i = 0; i < displace->GetNumberOfTuples();i++){
            pu_i = displace->GetTuple3(i);
            for (int j = 0; j < displace->GetNumberOfComponents();j++) {
                myfile << pu_i[j] << " ";
            }
            myfile << endl;
        }
        myfile.close();
    }
#endif


    return EXIT_SUCCESS;
}

