#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "TProfile.h"
#include <cstdlib>
#include <sstream>

using namespace std;

#define L 36814   //RunG : 7593


TFile *Target;
typedef vector<pair<TFile*, Double_t> > vec_pair;
typedef vector<pair<TFile*, Double_t> >::const_iterator vec_pair_it;

void MergeRootfile( TDirectory *target, const vector<pair<TFile*, Double_t> >& vFileList);

void haddws_yesWei_noRand() {
  /**********************************************************
   * in an interactive ROOT session, edit the file names,
  * corresponding weights, and target name. Then
   * root:> .x haddws.C+
  **********************************************************/

   Target = TFile::Open( "AccXEff_yesWei_noRand.root", "RECREATE" );

  vec_pair vFileList;
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy120.root"), (1.975*L)/2977.600 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy200.root"), (19.32*L)/100000 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy400.root"), (2.731*L)/100000 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy800.root"), (0.241*L)/98400 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy1400.root"), (0.01678*L)/100000 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy2300.root"), (0.00139*L)/95106 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy3500.root"), (0.00008948*L)/100000 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy4500.root"), (0.000004135*L)/100000 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dy6000.root"), (4.56e-7*L)/100000 ));
  vFileList.push_back(make_pair(TFile::Open("histos_yesWei_noRand_dyInf.root"), (2.066e-8*L)/100000 ));


  for (vec_pair_it it = vFileList.begin(); it != vFileList.end(); ++it)
  {
    cout << "File/weight = " << it->first->GetName() << "/" << it->second << endl;
  }

   MergeRootfile( Target, vFileList );
  Target->Close();

  for (vec_pair_it it = vFileList.begin(); it != vFileList.end(); ++it)
  {
    delete (it->first);
  }
}

void MergeRootfile( TDirectory *target, const vector<pair<TFile*, Double_t> >& vFileList) {

   //  cout << "Target path: " << target->GetPath() << endl;
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );

  vec_pair_it it = vFileList.begin();
   TFile *first_source = (*it).first;
   const Double_t first_weight = (*it).second;
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   // loop over all keys in this directory
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {

      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();

      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
         // descendant of TH1 -> merge it

         //      cout << "Merging histogram " << obj->GetName() << endl;
         TH1 *h1 = (TH1*)obj;  h1->Sumw2();
      if (first_weight>0) {
        h1->Scale(first_weight);
      }

         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
      for (vec_pair_it nextsrc = vFileList.begin()+1; nextsrc != vFileList.end(); ++nextsrc ) {

            // make sure we are at the correct directory level by cd'ing to path
            (*nextsrc).first->cd( path );
        const Double_t next_weight = (*nextsrc).second;
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
               TH1 *h2 = (TH1*)key2->ReadObj();  h2->Sumw2();
          if (next_weight>0) h2->Scale(next_weight);
               h1->Add( h2 );
               delete h2;
            }
         }

      } else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

         // loop over all source files create a chain of Trees "globChain"
         const char* obj_name= obj->GetName();

         globChain = new TChain(obj_name);
         globChain->Add(first_source->GetName());
      for (vec_pair_it nextsrc = vFileList.begin()+1; nextsrc != vFileList.end(); ++nextsrc ) {
            globChain->Add(nextsrc->first->GetName());
         }

      } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
         // it's a subdirectory

         cout << "Found subdirectory " << obj->GetName() << endl;

         // create a new subdir of same name and title in the target file
         target->cd();
         TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

         // newdir is now the starting point of another round of merging
         // newdir still knows its depth within the target file via
         // GetPath(), so we can still figure out where we are in the recursion
         MergeRootfile( newdir, vFileList);

      } else {

         // object is of no type that we know or can handle
         cout << "Unknown object type, name: "
         << obj->GetName() << " title: " << obj->GetTitle() << endl;
      }

      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      if ( obj ) {
         target->cd();

         //!!if the object is a tree, it is stored in globChain...
         if(obj->IsA()->InheritsFrom( TTree::Class() ))
            globChain->Merge(target->GetFile(),0,"keep");
         else
            obj->Write( key->GetName() );
      }

   } // while ( ( TKey *key = (TKey*)nextkey() ) )

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}

