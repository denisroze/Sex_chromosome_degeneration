#include "sexchr.h"
#include <cmath>
using namespace std;

extern MTRand rnd;

// male fitness:

double Wmale(ind &parent, double expdom, double smax, double Qopt, double I, int nbS)
{
    int j, allele1, allele2;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus,explvscaled,transmale;
    
    transmale = (((parent.transm[0] > 0) ? parent.transm[0] : 0) + ((parent.transm[1] > 0) ? parent.transm[1] : 0))/2;
    
    // multiplicative effects of the fitness loci:

    for (j = 0; j < nbS; j++)
    {
        allele1 = j;
        allele2 = nbS+j;
        
        e1 = ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);  // cis_y
        e2 = ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0) ;  // cis_x

        explv = (e1 + e2) * transmale; //  expression level
        
        if (explv == 0)
            w *= 1 - smax;  // fitness decreased by 1-smax if no expression of gene

        else
        {
            // stabilizing selection on expression level:
            
            explvscaled= log(explv/Qopt);
            effectLocus = 1-smax*(1-exp(-I*explvscaled*explvscaled));

            if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                    h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e1/(e1+e2)), expdom);

                // fitness = w1 + h*(w2-w1)
                //w += parent.gene[nbS+j] + h * (parent.gene[j] - parent.gene[nbS+j]); <-- additive fitness
                effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2]));
            }
            else // if first allele is fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                    h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e2/(e1+e2)), expdom);

                // fitness = w1 + h*(w2-w1)
                //w += parent.gene[j] + h * (parent.gene[nbS+j] - parent.gene[j]); <-- additive fitness
                effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1]));
            }
            if (effectLocus < 1 - smax)
                w *= 1 - smax;
            else
                w *= effectLocus;
        }
    }
    return w;
}

//female fitness:

double Wfemale(ind &parent, double expdom, double smax, double Qopt, double I, int nbS)
{
    int j;
    double w = 1.0;
    double e1, e2, c1, c2, explv, h, effectLocus,explvscaled,transfemale;
    int allele1, allele2;

    transfemale = (((parent.transf[0] > 0) ? parent.transf[0] : 0) + ((parent.transf[1] > 0) ? parent.transf[1] : 0))/2;

    for (j = 0; j < nbS; j++)
    {
        allele1 = j;
        allele2 = nbS+j;
            
        e1 = ((parent.cis[allele1] > 0) ? parent.cis[allele1] : 0);
        e2 = ((parent.cis[allele2] > 0) ? parent.cis[allele2] : 0);
        explv = (e1 + e2) * transfemale ;

        if (explv == 0)
            w *= 1 - smax;

        else
        {
            // stabilizing selection on expression level:
            explvscaled= log(explv/Qopt);

            effectLocus = 1-smax*(1-exp(-I*explvscaled*explvscaled));      // NEW


            if (parent.gene[allele2] >= parent.gene[allele1]) // if second allele fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e1/(e1+e2)), expdom);

                effectLocus *= (parent.gene[allele2] + h * (parent.gene[allele1] - parent.gene[allele2])); // product fitnesses
            }
            else // if first allele fittest
            {
                // h is the dominance coefficient, H0 = 0.25
                //                   h = pow((e1 / (e1 + e2)), (-log(H0)/log(2)));
                h = pow((e2/(e1+e2)), expdom);

                effectLocus *= (parent.gene[allele1] + h * (parent.gene[allele2] - parent.gene[allele1])); // product fitnesses
            }
            if (effectLocus < 1 - smax)
                w *= 1 - smax;
            else
                w *= effectLocus;
        }
    }
    
    return w;
}



