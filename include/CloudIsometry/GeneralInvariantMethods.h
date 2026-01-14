#ifndef GENERALINVARIANTMETHODS_H
#define GENERALINVARIANTMETHODS_H


class GeneralInvariantMethods
{
public:
    GeneralInvariantMethods() {}
    virtual ~GeneralInvariantMethods() {}

    static double permutationSign(std::vector<int> & a)
    {
        int n = a.size();
        int cnt=0;
        for(int i=0; i<n; i++)
            for(int j=i+1; j<n; j++)
                if (a[i]>a[j]) cnt++;

        if (cnt%2 == 0)
        {
            return 1;
        }
        else
        {
            return -1;

        }

    }



protected:

private:
};

#endif // GENERALINVARIANTMETHODS_H
