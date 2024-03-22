void getMeans(TTree* t, std::string parName){

    double par;
    double avg = 0.0;
    t->SetBranchAddress(parName.c_str(), &par);

    for (int i=0; i<t->GetEntries();i++){
        t->GetEntry(i);
        avg += par;
    }

    avg = avg / t->GetEntries();
    printf("avg %s = %f\n", parName.c_str(), par);

}
