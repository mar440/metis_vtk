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
// MPI
//#include "mpi.h"
//#include <openmpi.h>
#include <set>
#include <stdexcept>
#include "metis.h"



using namespace std;
bool asciiOrBinaryVtu = true;


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
    vtkSmartPointer<vtkXMLUnstructuredGridReader> imesh =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    imesh->SetFileName(filename.c_str());
    imesh->Update();

    int _i, nCellArr = imesh->GetNumberOfCellArrays();
    string str_i;
    for (_i = 0; _i < nCellArr ; _i++) {
        str_i = imesh->GetCellArrayName(_i);
        if (str_i.compare("PartitionId") == 0){
            break;
        }
    }
    cout << str_i << endl;
    int nb_node = imesh->GetNumberOfPoints();
    int nb_elmt = imesh->GetNumberOfCells();
    cout << "imesh->GetNumberOfPoints();" << nb_node << endl;
    cout << "imesh->GetNumberOfCells();" << nb_elmt << endl;
    vtkDataArray * PartitionId = imesh->GetOutput()->GetCellData()->GetArray("PartitionId");
    vtkCellArray * offset = imesh->GetOutput()->GetCells();

    int ds = PartitionId->GetDataSize();

    cout << "PartitionId.GetDataSize() = " << ds << endl;
    cout << "nb_elmt                   = " << nb_elmt << endl;


    cout << "offset->GetDataSize() = " << offset->GetSize()  << endl;


    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    int nPi;
    int * eptr = new int[nb_elmt + 1];
    eptr[0] = 0;
    for (int i = 0 ; i < nb_elmt; i++){
        imesh->GetOutput()->GetCell(i,cell);
        nPi = cell->GetNumberOfPoints();
        eptr[i+1] = eptr[i] + nPi;
    }
    int * eind = new int[eptr[nb_elmt]];
    int cnt = 0;
    int nNodOnEl=0;
    int CellDim;
    for (int i = 0 ; i < nb_elmt; i++){
        imesh->GetOutput()->GetCell(i,cell);
        nPi = cell->GetNumberOfPoints();
        for (int j = 0; j < nPi; j++){
            eind[cnt] = cell->GetPointId(j);
            cnt ++;
        }
        if (CellDim <cell->GetCellDimension()){
            CellDim  = cell->GetCellDimension();
        }
    }

    int ncommon = 0;
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

    int nparts = 1;
    int objval;


    if (_nparts > 0 ){
       nparts = _nparts;
    }

    int *epart = new int [nb_elmt];
    int *npart = new int [eptr[nb_elmt]];

    METIS_PartMeshDual(&nb_elmt,        // number of elements in the mesh       Y
                       &nb_node,        //                                      Y
                       eptr,            //                                      Y
                       eind,            //                                      Y
                       (int *)NULL,     // vwgt                                 Y
                       (int *)NULL,     // vsize                                Y
                       &ncommon,        //                                      Y
                       &nparts,         //                                      N
                       (real_t*)NULL,  // tpwgts                                Y
                       options,         //                                      Y
                       &objval,         //                                      Y
                       epart,           //                                      N
                       npart);          //                                      Y


    double tuple[] = {0};
    for (int i = 0 ; i < nb_elmt; i++){
        tuple[0] = epart[i];
        PartitionId->SetTuple(i,tuple);
    }

    delete [] epart;
    delete [] npart;
    delete [] eptr;
    delete [] eind;

    string fname = "dmpFls/modifiedFile_" + to_string(nparts) + "_subs.vtu";
    vtkSmartPointer<vtkUnstructuredGrid> ModelMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ModelMesh->ShallowCopy(imesh->GetOutput());

    printVTK(ModelMesh,fname);

    return EXIT_SUCCESS;
}

