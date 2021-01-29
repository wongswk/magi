//
// Created by Shihao Yang on 6/5/19.
//

#include <pybind11/buffer_info.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <armadillo>

#include <tgtdistr.h>
#include <hmc.h>
#include <classDefinition.h>
#include "magi_main_py.h"


namespace py = pybind11;

typedef arma::vec ArmaVector;
typedef std::shared_ptr< ArmaVector > ArmaVectorPointer;
typedef arma::mat ArmaMatrix;
typedef std::shared_ptr < ArmaMatrix > ArmaMatrixPointer;
typedef arma::cube ArmaCube;
typedef std::shared_ptr < ArmaCube > ArmaCubePointer;

// const static bool COPY_BUFFER = false;
// const static auto ownership = py::return_value_policy::reference;

const static bool COPY_BUFFER = true;
const static auto ownership = py::return_value_policy::take_ownership;

arma::vec createVector(py::array_t< double >& array)
{
    auto vector = array.mutable_unchecked< 1 >();
    return arma::vec(&vector(0), vector.size(), COPY_BUFFER, true);
}

arma::mat createMatrix(py::array_t< double >& array)
{
    auto matrix = array.mutable_unchecked< 2 >();
    return arma::mat(&matrix(0, 0), matrix.shape(1), matrix.shape(0), COPY_BUFFER, true);
}

arma::cube createCube(py::array_t< double >& array)
{
    auto cube = array.mutable_unchecked< 3 >();
    return arma::cube(&cube(0, 0, 0), cube.shape(2), cube.shape(1), cube.shape(0), COPY_BUFFER, true);
}


PYBIND11_MODULE(pymagi, macro)
{
    /*
     * ARMADILLO FUNCTIONS AND TYPES
     *
     *
     *
     *
     */
    py::class_< ArmaVector, ArmaVectorPointer >(macro, "ArmaVector", py::buffer_protocol())
        .def(py::init< arma::uword >())
        .def(py::init(&createVector), ownership)
        .def("size", &ArmaVector::size)
        .def("__getitem__", [](ArmaVector& vector, const arma::uword i) { return py::cast(vector[i]); })
        .def("__setitem__",
             [](ArmaVector& vector, const arma::uword i, py::object value) {
                 vector[i] = value.cast< double >();
                 return py::cast(vector[i]);
             })
        .def(
            "__iter__",
            [](const ArmaVector& vector) { return py::make_iterator(vector.begin(), vector.end()); },
            py::keep_alive< 0, 1 >())
        .def_buffer([](ArmaVector& vector) {
            return py::buffer_info(vector.memptr(),
                                   sizeof(double),
                                   py::format_descriptor< double >::format(),
                                   1,
                                   {vector.size()},
                                   {sizeof(double)});
        })
        .def("__repr__", [](const ArmaVector& vector) {
            std::stringstream ss;
            ss << "[ ";
            for (arma::uword i = 0; i < vector.size(); i++)
            {
                ss << vector[i];
                if (i + 1 < vector.size())
                {
                    ss << ", ";
                }
            }
            ss << "]";
            return ss.str();
        });

    py::class_< ArmaMatrix, ArmaMatrixPointer >(macro, "ArmaMatrix", py::buffer_protocol())
        .def(py::init< arma::uword, arma::uword >())
        .def(py::init(&createMatrix), ownership)
        .def_readonly("n_rows", &ArmaMatrix::n_rows)
        .def_readonly("n_cols", &ArmaMatrix::n_cols)
        .def(
            "col",
            [](const ArmaMatrix& matrix, arma::uword i) { return ArmaVector(matrix.col(i)); },
            py::return_value_policy::move)
        .def(
            "row",
            [](const ArmaMatrix& matrix, arma::uword i) { return ArmaVector(matrix.row(i).t()); },
            py::return_value_policy::move)
        .def(
            "t",
            [](const ArmaMatrix& matrix) { return ArmaMatrix(matrix.t()); },
            py::return_value_policy::move)
        .def("__getitem__",
             [](ArmaMatrix& matrix, std::tuple< arma::uword, arma::uword > coordinates) {
                 const auto& i = std::get< 0 >(coordinates);
                 const auto& j = std::get< 1 >(coordinates);
                 return py::cast(matrix(i, j));
             })
        .def("__setitem__",
             [](ArmaMatrix& matrix, std::tuple< arma::uword, arma::uword > coordinates, py::object value) {
                 const auto& i = std::get< 0 >(coordinates);
                 const auto& j = std::get< 1 >(coordinates);
                 matrix(i, j) = value.cast< double >();
                 return py::cast(matrix(i, j));
             })
        .def(
            "__iter__",
            [](ArmaMatrix& matrix) { return py::make_iterator(matrix.begin(), matrix.end()); },
            py::keep_alive< 0, 1 >())
        .def_buffer([](ArmaMatrix& matrix) {
            std::vector< size_t > shape = {matrix.n_rows, matrix.n_cols};
            std::vector< size_t > stride = {sizeof(double), sizeof(double) * matrix.n_rows};
            return py::buffer_info(
                matrix.memptr(), sizeof(double), py::format_descriptor< double >::format(), 2, shape, stride);
        })
        .def("__repr__", [](const ArmaMatrix& matrix) {
            std::stringstream ss;
            ss << "[ ";
            for (unsigned int i = 0; i < matrix.n_cols; i++)
            {
                ss << "[";
                for (unsigned int k = 0; k < matrix.n_rows; k++)
                {
                    ss << matrix(k, i);
                    if (k + 1 < matrix.n_rows)
                    {
                        ss << ", ";
                    }
                }
                ss << "]";
                if (i + 1 < matrix.n_cols)
                {
                    ss << ",";
                }
            }
            ss << "]";
            return ss.str();
        });

    py::class_< ArmaCube, ArmaCubePointer >(macro, "ArmaCube", py::buffer_protocol())
        .def(py::init< arma::uword, arma::uword, arma::uword >())
        .def(py::init(&createCube), ownership)
        .def_readonly("n_rows", &ArmaCube::n_rows)
        .def_readonly("n_cols", &ArmaCube::n_cols)
        .def_readonly("n_slices", &ArmaCube::n_slices)
        .def("slice", [](const ArmaCube& cube, arma::uword j) { return cube.slice(j); })
        .def("__getitem__",
             [](ArmaCube& cube, std::tuple< arma::uword, arma::uword, arma::uword > coordinates) {
                 const auto& i = std::get< 0 >(coordinates);
                 const auto& j = std::get< 1 >(coordinates);
                 const auto& k = std::get< 2 >(coordinates);
                 return py::cast(cube(i, j, k));
             })
        .def("__setitem__",
             [](ArmaCube& cube,
                std::tuple< arma::uword, arma::uword, arma::uword > coordinates,
                py::object value) {
                 const auto& i = std::get< 0 >(coordinates);
                 const auto& j = std::get< 1 >(coordinates);
                 const auto& k = std::get< 2 >(coordinates);
                 cube(i, j, k) = value.cast< double >();
                 return py::cast(cube(i, j, k));
             })
        .def(
            "__iter__",
            [](ArmaCube& cube) { return py::make_iterator(cube.begin(), cube.end()); },
            py::keep_alive< 0, 1 >())
        .def_buffer([](ArmaCube& cube) {
            std::vector< size_t > shape = {cube.n_rows, cube.n_cols, cube.n_slices};
            std::vector< size_t > stride = {
                sizeof(double), sizeof(double) * cube.n_rows, sizeof(double) * cube.n_cols * cube.n_rows};
            return py::buffer_info(
                cube.memptr(), sizeof(double), py::format_descriptor< double >::format(), 3, shape, stride);
        })
        .def("__repr__", [](const ArmaCube& cube) {
            std::stringstream ss;
            ss << "[";
            for (unsigned int i = 0; i < cube.n_slices; i++)
            {
                ss << "s :[";
                for (unsigned int j = 0; j < cube.n_cols; j++)
                {
                    ss << "c [";
                    for (unsigned int k = 0; k < cube.n_rows; k++)
                    {
                        ss << cube(k, j, i);
                        if (k + 1 < cube.n_rows)
                        {
                            ss << ", ";
                        }
                    }
                    ss << "]";
                    if (j + 1 < cube.n_cols)
                    {
                        ss << ",";
                    }
                }
                ss << "]";
                if (i + 1 < cube.n_slices)
                {
                    ss << ",";
                }
                ss << "\n";
            }
            ss << "]";
            return ss.str();
        });


    /*
     * cpp classes
     */
    py::class_< lp >(macro, "lp")
        .def(py::init<>())
        .def_readwrite("value", &lp::value)
        .def_readwrite("gradient", &lp::gradient);

    py::class_< hmcstate >(macro, "hmcstate")
        .def(py::init<>())
        .def_readwrite("final", &hmcstate::final)
        .def_readwrite("finalp", &hmcstate::finalp)
        .def_readwrite("step", &hmcstate::step)
        .def_readwrite("trajH", &hmcstate::trajH)
        .def_readwrite("lprvalue", &hmcstate::lprvalue)
        .def_readwrite("apr", &hmcstate::apr)
        .def_readwrite("delta", &hmcstate::delta)
        .def_readwrite("acc", &hmcstate::acc)
        .def_readwrite("trajq", &hmcstate::trajq)
        .def_readwrite("trajp", &hmcstate::trajp);

    /*
     * cpp functions
     */
    macro.def(
        "phisigllik",
        &phisigllik,
        "",
        py::arg("phisig"),
        py::arg("yobs"),
        py::arg("dist"),
        py::arg("kernel"));

    py::class_< gpcov >(macro, "gpcov")
            .def(py::init<>())
            .def_readwrite("C", &gpcov::C)
            .def_readwrite("Cinv", &gpcov::Cinv)
            .def_readwrite("mphi", &gpcov::mphi)
            .def_readwrite("Kphi", &gpcov::Kphi)
            .def_readwrite("Kinv", &gpcov::Kinv);

    /*
     * cpp class with functionals
     */
    py::class_< OdeSystem >(macro, "OdeSystem")
        .def(py::init<>())
        .def_readwrite("fOde", &OdeSystem::fOde)
        .def_readwrite("fOdeDx", &OdeSystem::fOdeDx)
        .def_readwrite("fOdeDtheta", &OdeSystem::fOdeDtheta)
        .def_readwrite("name", &OdeSystem::name)
        .def_readwrite("thetaLowerBound", &OdeSystem::thetaLowerBound)
        .def_readwrite("thetaUpperBound", &OdeSystem::thetaUpperBound)
        .def_readwrite("thetaSize", &OdeSystem::thetaSize)
        .def_readwrite("xLowerBound", &OdeSystem::xLowerBound)
        .def_readwrite("xUpperBound", &OdeSystem::xUpperBound);

    /*
     * cpp function with functional input
     * https://pybind11.readthedocs.io/en/stable/advanced/cast/functional.html
     */
    macro.def(
        "basic_hmcC",
        &basic_hmcC,
        "",
        py::arg("lpr"),
        py::arg("initial"),
        py::arg("step"),
        py::arg("lb"),
        py::arg("ub"),
        py::arg("nsteps"),
        py::arg("traj"));

    macro.def(
        "generalMaternCov",
        &generalMaternCov,
        "",
        py::arg("phi"),
        py::arg("distSigned"),
        py::arg("complexity"));

    py::class_< MagiSolver >(macro, "MagiSolver")
        .def_readwrite("phiAllDimensions", &MagiSolver::phiAllDimensions)
        .def_readwrite("sigmaInit", &MagiSolver::sigmaInit)
        .def_readwrite("xInit", &MagiSolver::xInit)
        .def_readwrite("thetaInit", &MagiSolver::thetaInit)
        .def_readwrite("covAllDimensions", &MagiSolver::covAllDimensions)
        .def_readwrite("llikxthetasigmaSamples", &MagiSolver::llikxthetasigmaSamples);

    macro.def(
        "solveMagiPy",
        &solveMagiPy,
        "",
        py::arg("yFull"),
        py::arg("odeModel"),
        py::arg("tvecFull"),
        py::arg("sigmaExogenous"),
        py::arg("phiExogenous"),
        py::arg("xInitExogenous"),
        py::arg("thetaInitExogenous"),
        py::arg("muExogenous"),
        py::arg("dotmuExogenous"),
        py::arg("priorTemperatureLevel"),
        py::arg("priorTemperatureDeriv"),
        py::arg("priorTemperatureObs"),
        py::arg("kernel"),
        py::arg("nstepsHmc"),
        py::arg("burninRatioHmc"),
        py::arg("niterHmc"),
        py::arg("stepSizeFactorHmc"),
        py::arg("nEpoch"),
        py::arg("bandSize"),
        py::arg("useFrequencyBasedPrior"),
        py::arg("useBand"),
        py::arg("useMean"),
        py::arg("useScalerSigma"),
        py::arg("useFixedSigma"),
        py::arg("verbose"));

}

