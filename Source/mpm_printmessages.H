void PrintWelcomeMessage()
{
        amrex::Print() << " ===============================================\n";
        amrex::Print() << "        Welcome to EXAGOOP MPM Solver           \n";
        amrex::Print() << "        Developed by HPACFers at NREL           \n";
        amrex::Print() << "                 -Hari, Sree and Marc           \n";
        amrex::Print() << " ===============================================\n";
}

void PrintMessage(std::string msg,int print_length,bool begin)
{
        if(begin==true)
        {
                msg.append(print_length - msg.length(), '-');
                msg.append(1, '>');
                amrex::Print() <<msg;
        }
        else
        {
                msg=" Done";
                amrex::Print() <<msg;
        }
}

void PrintMessage(std::string msg,int print_length,bool begin,char c)
{
        if(begin==true)
        {
                msg.append(print_length - msg.length(), c);
                amrex::Print() <<msg;
        }
        else
        {
                msg=" Done";
                amrex::Print() <<msg;
        }
}

void PrintEndMessage(std::string msg)
{
  amrex::Print() <<msg;
}
