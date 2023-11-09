#include "Template_Classes.hpp"

void DomainForceFieldUnit::enqueue_calculate_force_on_boundaries() { // calculate forces from fluid on TYPE_S nodes
	kernel_calculate_force_on_boundaries.set_parameters(2u, t).enqueue_run(); // опнакелю - опхмхлюер оюпюлерпш наыецн назейрю
}