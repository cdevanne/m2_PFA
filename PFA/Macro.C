void Macro(double gapSize) 
    {
        ofstream myfile;
        myfile.open ("datas.txt", ios::app);
        myfile << gapSize <<"   ";
        myfile.close();
        tree->Process("SDHCAL_PFA.C");
    }

