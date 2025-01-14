#ifndef _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_
#define _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_ 1

#include <thor_scsi/core/elements_basis.h>


namespace thor_scsi {
	namespace elements {
		/**

		   Empty space between two "typical accelerator components"
		 */
		using thor_scsi::core::ElemType;
		class DriftType : public ElemType {
		public:
			inline DriftType(const Config &config) : ElemType(config){
				// transformation done by transfrom
				// ... done by Elemtype initialisation
				// ... pleonamsmus
			}


			const char* type_name() const override final { return "Drift"; };
			/**
			 *
			 * Todo:
			 *     replace function with mv operator
			 */
			//virtual void assign(const ElementVoidBase *other) override{
			//	const DriftType *O = static_cast<const DriftType*>(other);
			//	// transform = O->transform
			//	ElementVoidBase::assign(other);
			//}
			// double GetPB(const int n) { return 0e0; };

			inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<double> &ps) override final
			{ _pass(conf, ps); };
			inline void pass(thor_scsi::core::ConfigType &conf, ss_vect<tps> &ps) override final
			{ _pass(conf, ps); };

		private:
			template<typename T>
			void _pass(const thor_scsi::core::ConfigType &conf, ss_vect<T> &ps);
		};
	}
}

#endif // _THOR_SCSI_CORE_ELEMENTS_DRIFT_H_
/*
 * Local Variables:
 * mode: c++
 * c++-file-style: "python"
 * End:
 */
