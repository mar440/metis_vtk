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
#include <vtkIdTypeArray.h>

using namespace std;
bool asciiOrBinaryVtu = true;


//#ifndef VTK_CREATE
// // To create a smart pointer
//#define VTK_CREATE(type, name) \
//    vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
//#endif


void printVTK(vtkSmartPointer<vtkUnstructuredGrid> _unstructuredGrid,string fname){
    cout <<"vtkVersion: " << vtkVersion::GetVTKMajorVersion() << "\n";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fname.c_str());
    if (asciiOrBinaryVtu){
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

   // mesh_local->AddArray(PartitionId);

    delete [] epart;
    delete [] npart;
    delete [] eptr;
    delete [] eind;

}


//vtkUnstructuredGrid* ThresoldSelection(vtkUnstructuredGrid* usg, std::string fieldName, std::vector<int> values)
//{
//    if (values.size() == 1)
//    {
//        //WriteUG(usg, "C:\\Users\\a0h72255\\Desktop\\test.vtu");
//        VTK_CREATE(vtkThreshold, threshold);
//        threshold->SetInputData(usg);
//        threshold->ThresholdBetween(values[0] - 0.1, values[0] + 0.1);
//        threshold->SetInputArrayToProcess(0, 0, 0, 1, fieldName.c_str());
//        threshold->Update();
//
//        vtkUnstructuredGrid* res = vtkUnstructuredGrid::New();
//        res->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(threshold->GetOutput()));
//        return res;
//    }
//    // else make an append filter of all the material ids
//    else
//    {
//        VTK_CREATE(vtkAppendFilter, aFilter);
//        aFilter->MergePointsOff();
//        for (auto i : values)
//        {
//            VTK_CREATE(vtkThreshold, threshold);
//            threshold->SetInputData(usg);
//            threshold->ThresholdBetween(i - 0.1, i + 0.1);
//            threshold->SetInputArrayToProcess(0, 0, 0, 1, fieldName.c_str());
//            threshold->Update();
//
//            aFilter->AddInputData(threshold->GetOutput());
//        }
//
//        aFilter->Update();
//        vtkUnstructuredGrid* res = vtkUnstructuredGrid::New();
//        res->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(aFilter->GetOutput()));
//        return res;
//    }
//}

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

    vtkDataArray * PartitionId;
    if (flagPartitionId){
         PartitionId = meshGlobal->GetOutput()->GetCellData()->GetArray("PartitionId");
    }
    else{
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
    connectivityFilter->SetInputConnection(meshGlobal->GetOutputPort());
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

    int nAverCellsPerRegion = nCellsGlobal / nRegions;


    for (int i = 0; i < nRegions; i ++) {
        vtkSmartPointer <vtkUnstructuredGrid> omega_i = vtkSmartPointer< vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
        threshold->SetInputData(connectivityFilter->GetOutput());
        threshold->ThresholdBetween((double)i, (double)i);
        threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
        threshold->Update();
        omega_i->ShallowCopy(threshold->GetOutput());

        nCellsPerRegion[i] = omega_i->GetNumberOfCells();
        nPartsPerRegion[i] = double( nCellsPerRegion[i]) / nAverCellsPerRegion;
    }

    for (int i = 0; i < nRegions; i++){
        cout <<"=========================================="<< endl;
        cout <<"nPartsPerRegion  = " << nPartsPerRegion[i] << endl;
        cout <<"=========================================="<< endl;
    }



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
        decomposeRegion(omega_i, PartitionId, _nparts, prevPartitionId);
        prevPartitionId += _nparts;

        //  string fname2 = "dmpFls/Omega_"+to_string(i)+".vtu";
        //  printVTK(omega_i,fname2);
    }



//#if 0

    string fname = "dmpFls/modifiedFile_" + to_string(_nparts) + "_subs.vtu";
    vtkSmartPointer<vtkUnstructuredGrid> ModelMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ModelMesh->ShallowCopy(meshGlobal->GetOutput());
    printVTK(ModelMesh,fname);
//#endif

    if (!flagPartitionId){
       PartitionId->Delete();
    }
    return EXIT_SUCCESS;
}

