#include "../../headers/Elements/Element.h"

/*
            Invariant-Based Spectral Constitutive Tensor

            This file implements the constitutive tensor obtained from the
            invariant-based spectral decomposition of the strain tensor,
            as presented in:

            Nardi, D.C., Ferreira, A.R., Leonel, E.D.,
            "Revisiting Miehe’s Spectral Split: A basis-independent energy decomposition model for phase-field fracture", 2026.

            The implementation is provided as a modular template that can be
            adapted to different finite element frameworks.

            Data entry for the method:
            - int dim: spatial dimension (2 or 3);
            - double _divU : divergence of the displacement field;
            - double _degradFunc : degradation function value at the current iteration;
            - double _gradU[dim][dim] : gradient of the displacement field at the current iteration;
            - double D[dim][dim][dim][dim] : output constitutive tensor.

            Output:
            - double D[dim][dim][dim][dim] : constitutive tensor computed from the invariant-based spectral decomposition of the strain tensor.

            Any doubts or suggestions, please contact me at deborahnardi@usp.br
*/

template <int dim>
void Solid<dim>::getSpectralConstitutiveTensor(double _divU, double _degradFunc, double _gradU[dim][dim], double D[dim][dim][dim][dim])
{
    // ====================================================================================================================
    // ================================================= INITIALIZE =======================================================
    const PetscReal lame{material->getLameConstant()}; // get the Lame constant from the material properties
    const PetscReal mu{material->getShearModulus()};   // get the shear modulus from the material properties
    const PetscReal poisson{material->getPoisson()};   // get the Poisson's ratio from the material properties
    std::string caseName;                              // Variable to store the case name for post-processing purposes

    PetscReal strain[3][3]{};
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            strain[i][j] = 0.5 * (_gradU[i][j] + _gradU[j][i]);
    if (material->getPlaneAnalysis() == PLANE_STRESS)
        strain[2][2] = -(poisson / (1.0 - poisson)) * (strain[0][0] + strain[1][1]);
    // ====================================================================================================================
    // ====================================================================================================================
    const PetscReal volStrain{(1.0 / 3.0) * trace<3>(strain)};
    PetscReal deviatoricStrain[3][3]{};
    getDeviatoric<3>(strain, deviatoricStrain); // Get the deviatoric strain tensor
    const PetscReal eeq{sqrt(2.0 / 3.0 * doubleContraction<3>(deviatoricStrain, deviatoricStrain))};

    //  =====================================================================================================================
    //  ============================================== NUMERICAL NORMALIZATION ==============================================
    //  =====================================================================================================================

    /*
        This numerical normalization is not strictly necessary, but it significantly reduces the number of repetead operations,
        as the normalized strain tensor and its invariants are used multiple times in the computation of the constitutive tensor.
        Basically, the strain tensor is normalized by the equivalent strain, and the constitutive tensor is computed in terms of
        the normalized strain tensor (see the operation named "product" below).
    */
    PetscReal a[3][3]{};
    product<3>(strain, 1.0 / eeq, a);             // Normalize the strain tensor
    const PetscReal volStrainA{trace<3>(a) / 3.}; // Normalized volumetric strain
    PetscReal devA[3][3]{};                       // Normalized deviatoric strain tensor
    getDeviatoric<3>(a, devA);                    // Get the deviatoric normalized strain tensor
    PetscReal devAQuad[3][3]{};                   // Normalized Quadratic deviatoric strain tensor
    product<3>(devA, devA, devAQuad);             // Get the normalized quadratic deviatoric strain tensor

    const PetscReal detddA{getMatrixDeterminant<3>(devA)};
    // =====================================================================================================================
    // =====================================================================================================================
    PetscReal etaReal{0.}, cos3eta{0.};
    if (std::abs(eeq) > 1.e-12) // or 1.e-10
    {
        cos3eta = 4. * detddA;
        cos3eta = std::clamp(cos3eta, -1.0, 1.0);
        etaReal = std::acos(cos3eta) / 3.0;
    }
    const PetscReal etaStar[]{etaReal,
                              etaReal - 2.0 * M_PI / 3.0,
                              etaReal + 2.0 * M_PI / 3.0}; // Lode angle
    const PetscReal principalValues[]{volStrain + eeq * std::cos(etaStar[0]),
                                      volStrain + eeq * std::cos(etaStar[1]),
                                      volStrain + eeq * std::cos(etaStar[2])};

    const PetscReal pTol{1.e-8};
    bool strainIsNull{(eeq < pTol)};
    //  =====================================================================================================================
    //  ========================================== INITIALIZE DERIVATIVES ===================================================
    //  =====================================================================================================================
    PetscReal d_ep_de[3][3][3]{}, d2_ep_de2[3][3][3][3][3]{};            // First and second derivative of the energy functional with respect to the strain tensor
    PetscReal DPlus[dim][dim][dim][dim]{}, DMinus[dim][dim][dim][dim]{}; // Constitutive tensors for the positive and negative principal values
    // d_eeq_deij--------------------------------------------
    PetscReal d_eeq_de[3][3]{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            d_eeq_de[i][j] = 2.0 / 3.0 * devA[i][j];

    //  d2_eeq_deij_dekl -------------------------------------------------
    PetscReal d2_eeq_de2[3][3][3][3]{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            d2_eeq_de2[i][j][i][j] += 1. / 3.;
            d2_eeq_de2[i][j][j][i] += 1. / 3.;
            d2_eeq_de2[i][i][j][j] -= 2. / 9.;
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    d2_eeq_de2[i][j][k][l] -= 4. / 9. * devA[i][j] * devA[k][l];

    //  vij ---------------------------------------------
    PetscReal v[3][3]{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            v[i][j] = 2. * devAQuad[i][j] - 4. * detddA * devA[i][j];
    for (int i = 0; i < 3; i++)
        v[i][i] -= 1.0;

    //  dv_deij---------------------------------------------------
    PetscReal dv_de[3][3][3][3]{};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                dv_de[i][j][i][k] += devA[j][k];
                dv_de[i][j][k][i] += devA[k][j];
                dv_de[i][j][k][j] += devA[i][k];
                dv_de[i][k][k][j] += devA[i][j];

                dv_de[i][j][k][k] += (2. - 4. / 3.) * devA[i][j];

                dv_de[i][i][k][j] -= 8. / 3. * devA[k][j];
            }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            dv_de[i][j][i][j] -= 2. * detddA;
            dv_de[i][j][j][i] -= 2. * detddA;
            dv_de[i][i][j][j] += 4. / 3. * detddA;
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                {
                    dv_de[i][j][k][l] += 8. / 3. * devAQuad[i][j] * devA[k][l];
                    dv_de[i][j][k][l] -= 4. * devAQuad[k][l] * devA[i][j];
                }
    //-------------------------------------------------------------
    PetscReal tol{1.e-12}; //-8
    if (strainIsNull)
    {
        caseName = "NULL";
        const PetscReal c = 2.0 * mu;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                D[i][i][j][j] += lame;
                D[i][j][i][j] += _degradFunc * c;
            }
    }
    else if (std::fabs(principalValues[1] - principalValues[0]) > tol &&
             std::fabs(principalValues[2] - principalValues[1]) > tol)
    {
        // =====================================================================================================================
        //                                       GENERAL CASE: 3 DISTINCT PRINCIPAL VALUES
        // =====================================================================================================================
        caseName = "GENERAL";
        const PetscReal sin3eta{std::sin(3.0 * etaStar[0])};
        const PetscReal cos3eta{std::cos(3.0 * etaStar[0])};

        PetscReal d_eta_de[3][3]{};
        const PetscReal aux2{-2. / (3. * sin3eta)};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                d_eta_de[i][j] = aux2 * (2. * devAQuad[i][j] - 4. * detddA * devA[i][j]);
        for (int i = 0; i < 3; i++)
            d_eta_de[i][i] -= aux2;

        const PetscReal aux3{2.0 / 3.0 * (1.0 / sin3eta)};
        const PetscReal aux4{-9.0 / 2.0 * cos3eta};
        PetscReal d2_eta_de2[3][3][3][3]{};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        d2_eta_de2[i][j][k][l] = aux3 * (aux4 * d_eta_de[i][j] * d_eta_de[k][l] + 10. / 3. * devA[k][l] * v[i][j] - dv_de[i][j][k][l]);
        //  =====================================================================================================================
        // COMPUTE FIRST AND SECOND DERIVATIVES
        for (int ipp = 0; ipp < 3; ipp++)
        {
            const PetscReal cosEta{std::cos(etaStar[ipp])};
            const PetscReal sinEta{std::sin(etaStar[ipp])};

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    d_ep_de[ipp][i][j] += d_eeq_de[i][j] * cosEta - sinEta * d_eta_de[i][j];
            for (int i = 0; i < 3; i++)
                d_ep_de[ipp][i][i] += 1. / 3.;

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            d2_ep_de2[ipp][i][j][k][l] = 1. / eeq *
                                                         ((d2_eeq_de2[i][j][k][l] - d_eta_de[i][j] * d_eta_de[k][l]) * cosEta -
                                                          (d_eeq_de[i][j] * d_eta_de[k][l] + d_eta_de[i][j] * d_eeq_de[k][l] + d2_eta_de2[i][j][k][l]) * sinEta);
        }
        //  =====================================================================================================================
        // COMPUTE THE CONSTITUTIVE TENSOR
        const PetscReal c{2.0 * mu};
        for (int ipp = 0; ipp < 3; ipp++)
        {
            PetscReal(*Daux)[dim][dim][dim][dim] = (principalValues[ipp] > 0.0) ? &DPlus : &DMinus;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    for (int k = 0; k < dim; k++)
                        for (int l = 0; l < dim; l++)
                            (*Daux)[i][j][k][l] += c * (d_ep_de[ipp][i][j] * d_ep_de[ipp][k][l] + principalValues[ipp] * d2_ep_de2[ipp][i][j][k][l]);
        }
    }
    else
    {
        // ====================================================================================================================
        //                                  CHECK IF DELTA SHOULD'VE BE IMPOSED IN DIRECTION 22
        // ===================================================================================================================
        const PetscReal detEps{getMatrixDeterminant<3>(strain)};
        int alpha{0};

        if (strain[1][1] == strain[2][2] && detEps == strain[0][0] * strain[1][1] * strain[2][2] && strainIsNull == false)
            alpha{1};

        // =====================================================================================================================
        //                                      PARTICULAR CASE: 2 PRINCIPAL VALUES ARE EQUAL
        // =====================================================================================================================
        const PetscReal denum{sqrt(devA[alpha][alpha] * devA[alpha][alpha] + 2.0 * devAQuad[alpha][alpha] - 3.0 * devAQuad[alpha][alpha] * devAQuad[alpha][alpha])};
        const PetscReal coef{-1.0 / (3.0 * sqrt(3))};
        const PetscReal cc{coef / denum};

        PetscReal d_eta_de[3][3]{};
        PetscReal matAux[3][3]{};
        matAux[0][0] = 1.0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                d_eta_de[i][j] = cc *
                                 (4. * devAQuad[i][j] * devA[alpha][alpha] - 6. * devAQuad[alpha][alpha] * devA[i][j] +
                                  3. * (devA[i][alpha] * matAux[j][alpha] + devA[j][alpha] * matAux[i][alpha]) - 4. * devA[alpha][alpha] * matAux[i][j] + devA[i][j] -
                                  2. * detddA * (3. * matAux[i][alpha] * matAux[j][alpha] - matAux[i][j]));

        //  COMPUTE FIRST DERIVATIVE
        for (int ipp = 0; ipp < 3; ipp++)
        {
            const PetscReal cosEta{std::cos(etaStar[ipp])};
            const PetscReal sinEta{std::sin(etaStar[ipp])};

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    d_ep_de[ipp][i][j] += d_eeq_de[i][j] * cosEta - sinEta * d_eta_de[i][j];
            for (int i = 0; i < 3; i++)
                d_ep_de[ipp][i][i] += 1. / 3.;
        }
        // =====================================================================================================================
        if (abs(etaStar[0]) < tol)
        {
            caseName = "PARTICULAR I";
            /*
                cos(3eta) = 1; eta = 0
                IF eta = 0 -> ep2 = ep3
                We must check either if they are negative or positive;
                We must check the signal of ep1 as well.
            */
            PetscReal sum[3][3][3][3]{}, sumEp1[3][3][3][3]{};
            const PetscReal c2{1.0 / eeq}, c3{principalValues[0] / eeq};
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            sum[i][j][k][l] = c2 *
                                              ((d_eta_de[i][j] * d_eta_de[k][l] - d2_eeq_de2[i][j][k][l]) * principalValues[1] -
                                               (d_eta_de[i][j] * d_eta_de[k][l] + 2. / 9. * dv_de[i][j][k][l]) * principalValues[0]);
                            sumEp1[i][j][k][l] = c3 * (d2_eeq_de2[i][j][k][l] + 2.0 / 9.0 * dv_de[i][j][k][l]);
                        }
            //  =====================================================================================================================
            const PetscReal c{2.0 * mu};

            PetscReal(*Daux)[dim][dim][dim][dim] = (principalValues[1] > pTol && principalValues[2] > pTol) ? &DPlus : &DMinus;
            PetscReal(*Daux2)[dim][dim][dim][dim] = (principalValues[0] > pTol) ? &DPlus : &DMinus;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    for (int k = 0; k < dim; k++)
                        for (int l = 0; l < dim; l++)
                        {
                            (*Daux)[i][j][k][l] += c * (d_ep_de[1][i][j] * d_ep_de[1][k][l] + d_ep_de[2][i][j] * d_ep_de[2][k][l] + sum[i][j][k][l]);
                            (*Daux2)[i][j][k][l] += c * (d_ep_de[0][i][j] * d_ep_de[0][k][l] + sumEp1[i][j][k][l]);
                        }
        }
        else
        {
            caseName = "PARTICULAR II";
            /*
                eta = n Pi/3., where n is an integer (multiple of 60 degrees)
                IF eta = pi/3 -> ep1 = ep2
                    We must check either if they are negative or positive;
                    We must check the signal of ep3 as well.
            */
            PetscReal sum[3][3][3][3]{}, sumEp3[3][3][3][3]{};
            const PetscReal c2{1.0 / eeq}, c3{principalValues[2] / eeq};
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            sum[i][j][k][l] = c2 *
                                              ((d2_eeq_de2[i][j][k][l] - d_eta_de[i][j] * d_eta_de[k][l]) * principalValues[0] -
                                               (d_eta_de[i][j] * d_eta_de[k][l] - 2. / 9. * dv_de[i][j][k][l]) * principalValues[2]);
                            sumEp3[i][j][k][l] = c3 * (2.0 / 9.0 * dv_de[i][j][k][l] - d2_eeq_de2[i][j][k][l]);
                        }
            //  =====================================================================================================================
            const PetscReal c{2.0 * mu};
            PetscReal(*Daux)[dim][dim][dim][dim] = (principalValues[0] > pTol && principalValues[1] > pTol) ? &DPlus : &DMinus;
            PetscReal(*Daux2)[dim][dim][dim][dim] = (principalValues[2] > pTol) ? &DPlus : &DMinus;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    for (int k = 0; k < dim; k++)
                        for (int l = 0; l < dim; l++)
                        {
                            (*Daux)[i][j][k][l] += c * (d_ep_de[0][i][j] * d_ep_de[0][k][l] + d_ep_de[1][i][j] * d_ep_de[1][k][l] + sum[i][j][k][l]);
                            (*Daux2)[i][j][k][l] += c * (d_ep_de[2][i][j] * d_ep_de[2][k][l] + sumEp3[i][j][k][l]);
                        }
        }
    }

    if (!strainIsNull)
    {
        PetscReal(*Daux)[dim][dim][dim][dim] = (_divU > 0.0) ? &DPlus : &DMinus;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                (*Daux)[i][i][j][j] += lame;

        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                for (int k = 0; k < dim; k++)
                    for (int l = 0; l < dim; l++)
                        D[i][j][k][l] = _degradFunc * DPlus[i][j][k][l] + DMinus[i][j][k][l];
    }
}

template class Solid<2>;
template class Solid<3>;