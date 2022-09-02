void Macro(double gapSize) 
    {
        ofstream myfile;
        myfile.open ("datas.txt", ios::app);
        myfile << gapSize <<"   ";
        myfile.close();
        tree->Process("PFA_SDHCAL_TIMING.C");
    }

