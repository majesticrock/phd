#pragma once
#include "GlobalDefinitions.hpp"
#include <SymbolicOperators/TermLoader.hpp>
#include <SymbolicOperators/WickTerm.hpp>
#include <Utility/better_to_string.hpp>
#include "SCModel.hpp"
#include <memory>
#include <map>

#ifndef _complex
#define _XP
#endif

#ifndef _XP
#include <Utility/Numerics/iEoM/GeneralResolvent.hpp>
#define __ieom_algorithm Utility::Numerics::iEoM::GeneralResolvent<ModeHelper, c_complex>
#else
#include <Utility/Numerics/iEoM/XPResolvent.hpp>
#define __ieom_algorithm Utility::Numerics::iEoM::XPResolvent<ModeHelper, c_float>
#endif

namespace Continuum {
	class ModeHelper : public __ieom_algorithm
	{
	private:
		friend struct __ieom_algorithm;
		using _parent = __ieom_algorithm;
		using m_iterator = InnerIterator;

		c_float compute_momentum(SymbolicOperators::Momentum const& momentum, c_float k, c_float l, c_float q = 0) const;
		c_complex get_expectation_value(SymbolicOperators::WickOperator const& op, c_float momentum) const;

		c_complex compute_phonon_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const;
		c_complex compute_em_sum(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const;
	protected:
		SymbolicOperators::TermLoader wicks;
		//size_t TOTAL_BASIS{};
		constexpr static int hermitian_size = 2;
		constexpr static int antihermitian_size = 1;
		constexpr static int number_of_basis_terms = hermitian_size + antihermitian_size;

		static int hermitian_discretization;
		static int antihermitian_discretization;
		static int total_matrix_size;

		std::unique_ptr<SCModel> model;

		void createStartingStates();
		void fillMatrices();
		void fill_M();

		void fill_block_M(int i, int j);
		void fill_block_N(int i, int j);

		c_complex computeTerm(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const;
	public:
		std::vector<c_float> continuum_boundaries() const;

		SCModel& getModel() {
			return *model;
		};
		const SCModel& getModel() const {
			return *model;
		};
		ModeHelper(ModelInitializer const& init);
	};
}

#undef __ieom_algorithm