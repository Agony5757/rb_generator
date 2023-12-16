#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-value"
#endif

#include "rb22.h"
#include "rb44.h"

PYBIND11_MODULE(RBGeneratorCpp, m)
{
	py::enum_<Generator22Enum>(m, "Generator22Enum")
		.value("Generator_I22", Generator_I22)
		.value("Generator_X", Generator_X)
		.value("Generator_Y", Generator_Y)
		.value("Generator_SX", Generator_SX)
		.value("Generator_SY", Generator_SY)
		.value("Generator_SXdag", Generator_SXdag)
		.value("Generator_SYdag", Generator_SYdag)
		.export_values();

	py::class_<KeyValueClifford22>(m, "KeyValueClifford22")
		.def_readonly("arr", &KeyValueClifford22::arr)
		.def_readonly("buffer", &KeyValueClifford22::buffer)
		.def_readonly("size", &KeyValueClifford22::size)
		;	

	py::class_<RB22>(m, "RB22")
		.def(py::init<>())
		.def_readonly("N", &RB22::N)
		.def_readonly("pos_of_identity", &RB22::pos_of_identity)
		.def("load_from_file", &RB22::load_from_file)
		.def("get_multiplication_table_elem", &RB22::get_multiplication_table_elem)
		.def("get_inverse", &RB22::get_inverse)
		.def("get_group_elements", &RB22::get_group_elements, py::return_value_policy::reference)
		.def("get_table", &RB22::get_table, py::return_value_policy::reference)
		.def("get_inverse_table", &RB22::get_inverse_table, py::return_value_policy::reference)
		.def("get_matrix", &RB22::get_matrix, py::return_value_policy::reference)
		.def("get_generator", &RB22::get_generator, py::return_value_policy::reference)
		.def("get_generator_size", &RB22::get_generator_size)
		.def("get_full_sequence_and_inverse_sequence", &RB22::get_full_sequence_and_inverse_sequence)
		;

	py::class_<RB44>(m, "RB44")
		.def(py::init<>())
		.def_readonly("N", &RB44::N)
		.def_readonly("pos_of_identity", &RB44::pos_of_identity)
		.def("load_from_file", &RB44::load_from_file)
		.def("get_multiplication_table_elem", &RB44::get_multiplication_table_elem)
		.def("get_inverse", &RB44::get_inverse)
		.def("get_group_elements", &RB44::get_group_elements, py::return_value_policy::reference)
		.def("get_table", &RB44::get_table, py::return_value_policy::reference)
		.def("get_inverse_table", &RB44::get_inverse_table, py::return_value_policy::reference)
		.def("get_matrix", &RB44::get_matrix, py::return_value_policy::reference)
		.def("get_generator", &RB44::get_generator, py::return_value_policy::reference)
		.def("get_generator_size", &RB44::get_generator_size)
		.def("get_full_sequence_and_inverse_sequence", &RB44::get_full_sequence_and_inverse_sequence)
		;

	m.doc() = "[Module RBGeneratorCpp]";
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif