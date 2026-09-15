// Host-only formal evaluator handshake. Does not read/write populations or GPU metadata.
#include <chrono>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <unistd.h>
inline void formalExternalEvaluate(std::string const& run,int step) {
    const std::string tmp=run+"/diagnostics/evaluation_request.tmp";
    const std::string request=run+"/diagnostics/evaluation_request.txt";
    { std::ofstream out(tmp);out<<step<<'\n';out.flush();if(!out)throw std::runtime_error("RUNTIME_ERROR: evaluator request write"); }
    if(std::rename(tmp.c_str(),request.c_str())!=0)throw std::runtime_error("RUNTIME_ERROR: evaluator request rename");
    const auto begin=std::chrono::steady_clock::now();
    while(true) {
        {std::ifstream f(run+"/diagnostics/evaluation_error.txt");int failed=-1;f>>failed;if(f&&failed==step)throw std::runtime_error("RUNTIME_ERROR: external evaluator failed");}
        {std::ifstream f(run+"/diagnostics/decision.txt");int evaluated=-1,stop=-1,candidate=-1;f>>evaluated>>stop>>candidate;
            if(f&&evaluated==step){if((stop!=0&&stop!=1)||(candidate!=0&&candidate!=1))throw std::runtime_error("RUNTIME_ERROR: invalid evaluator decision");return;}
            if(f&&evaluated>step)throw std::runtime_error("RUNTIME_ERROR: evaluator future decision");
        }
        if(std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()>300.)
            throw std::runtime_error("RUNTIME_ERROR: evaluator service unavailable for300seconds");
        usleep(100000);
    }
}
