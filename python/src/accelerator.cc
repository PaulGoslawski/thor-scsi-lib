#include <pybind11/stl.h>
// #include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include "thor_scsi.h"
#include <thor_scsi/core/machine.h>
#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/std_machine/accelerator.h>

//namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace ts = thor_scsi;
namespace py = pybind11;


void py_thor_scsi_init_accelerator(py::module &m)
{

	// m.def("register_elements", &register_elements);
	// needs to be done only once
	ts::register_elements();

	py::class_<ts::Accelerator, std::shared_ptr<ts::Accelerator>>(m, "Accelerator")
		.def("find",                 &ts::Accelerator::find)
		// unique_ptr not working ... check for memory management ...
		.def("elementsWithName",     &ts::Accelerator::elementsWithName)
		// unique_ptr not working ... check for memory management ...
		.def("elementsWithNameType", &ts::Accelerator::elementsWithNameType)
		.def("__len__",              &ts::Accelerator::size)
		.def("__getitem__", py::overload_cast<size_t>(&ts::Accelerator::at))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_dbl&, size_t, int>(&ts::Accelerator::propagate))
		.def("propagate", py::overload_cast<tsc::ConfigType&, ts::ss_vect_tps&, size_t, int>(&ts::Accelerator::propagate))

		.def(py::init<const Config &>());



}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
