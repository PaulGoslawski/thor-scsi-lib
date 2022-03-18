#include <thor_scsi/std_machine/std_machine.h>
#include <thor_scsi/elements/drift.h>
#include <thor_scsi/elements/quadrupole.h>
#include <thor_scsi/elements/sextupole.h>
#include <thor_scsi/elements/octupole.h>
#include <thor_scsi/elements/marker.h>
#include <thor_scsi/elements/bpm.h>
#include <thor_scsi/elements/cavity.h>
#include <thor_scsi/elements/bending.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

int
register_elements(void)
{
	tsc::Machine::registerElement<tse::DriftType>("Drift");
	tsc::Machine::registerElement<tse::CavityType>("Cavity");
	tsc::Machine::registerElement<tse::MarkerType>("Marker");
	tsc::Machine::registerElement<tse::BPMType>("BPM");
	tsc::Machine::registerElement<tse::QuadrupoleType>("Quadrupole");
	tsc::Machine::registerElement<tse::SextupoleType>("Sextupole");
	tsc::Machine::registerElement<tse::OctupoleType>("Octupole");
	tsc::Machine::registerElement<tse::BendingType>("Bending");

	// tsc::Machine::registerElement<tse::MpoleType>("mpole");
	return 1;
}

/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
