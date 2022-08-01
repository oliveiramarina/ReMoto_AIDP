! '''
!     Neuromuscular simulator in Python.
!     Copyright (C) 2018  Renato Naville Watanabe

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

!     Contact: renato.watanabe@ufabc.edu.br
! '''

module SynapsesFactoryModule
    ! '''
    ! Module to build all the synapses in the system.
    ! '''
    use ConfigurationClass
    use MotorUnitPoolClass
    use InterneuronPoolClass
    use NeuralTractClass
    use DynamicalArrays
    use SynapsePointerClass
    use SynapticNoiseClass
    use CharacterMatrixClass
    use AfferentPoolClass
    implicit none
    integer, parameter :: wp = kind(1.0d0)
    private :: wp
    contains

        function synapseFactory(conf, neuralTractPools, &
            motorUnitPools, interneuronPools, afferentPools) result(synapticNoisePools)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pools**: list of all the pools in the system.

            
            ! '''
            class(Configuration), intent(inout) :: conf
            class(NeuralTract), intent(inout) :: neuralTractPools(:)
            class(AfferentPool), intent(inout) :: afferentPools(:)
            class(MotorUnitPool), intent(inout), target:: motorUnitPools(:)
            class(InterneuronPool), intent(inout), target:: interneuronPools(:)
            type(SynapticNoise), allocatable:: synapticNoisePools(:)
            type(SynapticNoise), allocatable:: tempSynNoise(:)
            integer :: numberOfSynapses, numberOfSynapticNoise
            integer :: poolOut, unitOut, synapseIn, poolIn, unitIn, compartmentIn, synapseComp
            character(len = 80) :: neuralSource, paramTag, paramChar
            real(wp) :: tau, var, gmax, conn, delay, declineFactor
            character(len=80) :: dyn
            integer :: pos, i, newIndex, j
            real(wp) :: randomNumber
            type(CharacterMatrix) :: NoiseSynapsesOut
            real(wp) :: neuronsDistance, weight
            type(SynapsePointer), dimension(:), allocatable:: tempTransmitSpikes
            

            ! ## Total number of synapses in the system.
            numberOfSynapses = 0
            
            !NeuralTract to MotorUnitPool
            if (size(neuralTractPools)>0 .and. size(motorUnitPools)>0) then
                do poolOut = 1, size(neuralTractPools)
                    do unitOut = 1, size(neuralTractPools(poolOut)%unit)
                        neuralSource = trim(neuralTractPools(poolOut)%pool) // '-'&
                         // neuralTractPools(poolOut)%unit(unitOut)%neuronKind
                        neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut = conf%determineSynapses(neuralSource)
                        do synapseIn = 1, size(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item)
                            paramTag = 'Con:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                            read(paramChar, *)conn
                            conn = conn / 100.0

                            paramTag = 'gmax:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)                            
                            paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                            read(paramChar,*)gmax

                            paramTag = 'delay:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, neuralTractPools(poolOut)%pool, 0)
                            read(paramChar, *)delay
                            
                            paramTag = 'dec:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, neuralTractPools(poolOut)%pool, 0)
                            if (trim(paramChar)=='inf') then
                                declineFactor = 1e6
                            else 
                                read(paramChar, *)declineFactor
                            end if
                            
                            paramTag = 'dyn:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            dyn = conf%parameterSet(paramTag, neuralTractPools(poolOut)%pool, 0)
                            
                            if (trim(dyn).ne.'None' .and. trim(dyn).ne.'none') then
                                paramTag = 'var:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                
                                paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                                read(paramChar, *)var

                                paramTag = 'tau:' // trim(neuralTractPools(poolOut)%pool) // '-' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,neuralTractPools(poolOut)%pool, 0)
                                read(paramChar, *)tau
                            else
                                var = 0.0
                                tau = 1e6
                            end if
                            
                            do poolIn = 1, size(motorUnitPools)
                                if (trim(neuralTractPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)==&
                                    trim(motorUnitPools(poolIn)%pool)) then
                                    do unitIn = 1, size(motorUnitPools(poolIn)%unit)
                                        do compartmentIn = 1, size(motorUnitPools(poolIn)%unit(unitIn)%Compartments)
                                            if ((trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(1)%string)==&
                                                trim(motorUnitPools(poolIn)%pool)) .and.&
                                                (trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(2)%string)==&
                                                trim(motorUnitPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                                (trim(neuralTractPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(3)%string)==&
                                                trim(motorUnitPools(poolIn)%unit(unitIn)%&
                                                Compartments(compartmentIn)%compKind))) then
                                                call random_number(randomNumber)
                                                if (randomNumber <= conn) then
                                                    do synapseComp = 1, size(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                            Compartments(compartmentIn)%SynapsesIn)
                                                        if (motorUnitPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%synapseKind ==&
                                                            neuralTractPools(poolOut)%unit(unitOut)%&
                                                            SynapsesOut%item(synapseIn)%item(4)%string) then
                                                            
                                                            if (declineFactor<1e5) then
                                                                neuronsDistance = abs(motorUnitPools(poolIn)%&
                                                                                unit(unitIn)%position_mm &
                                                                            - neuralTractPools(poolOut)%&
                                                                            unit(unitOut)%position_mm)
                                                                weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                                gmax = gmax * weight
                                                            end if
                                                            !gmax1 = gmax + unitIn*0.1/size(motorUnitPools(poolIn)%unit)
                                                            call motorUnitPools(poolIn)%unit(unitIn)%&
                                                                 Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                                 addConductance(gmax, delay, dyn, var, tau)
                                                            
                                                            if (allocated(neuralTractPools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)) then
                                                                
                                                                
                                                                allocate(tempTransmitSpikes(size(neuralTractPools(poolOut)%&
                                                                         unit(unitOut)%transmitSpikesThroughSynapses)))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call tempTransmitSpikes(j)%&
                                                                        assignSynapse(neuralTractPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                                end do
                                                                
                                                                
                                                                deallocate(neuralTractPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses)
                                                                
                                                                allocate(neuralTractPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call neuralTractPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(j)%&
                                                                    assignSynapse(tempTransmitSpikes(j)%synapse)
                                                                end do

                                                                neuralTractPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                        SynapsePointer()        

                                                                call neuralTractPools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                                assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                Compartments(compartmentIn)%SynapsesIn(synapseComp))  
                                                                
                                                                deallocate(tempTransmitSpikes)
                                                            else 
                                                                allocate(neuralTractPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1))
                                                                neuralTractPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)=&
                                                                        SynapsePointer()        

                                                                call neuralTractPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)%&
                                                                        assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                        Compartments(compartmentIn)%&
                                                                        SynapsesIn(synapseComp))
                                                            end if            
                                                            
                                                            newIndex = size(motorUnitPools(poolIn)%&
                                                                unit(unitIn)%Compartments(compartmentIn)%&
                                                                SynapsesIn(synapseComp)%gmax_muS)
                                                            
                                                            call integerAddToList(neuralTractPools(poolOut)%&
                                                            unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                            
                                                            numberOfSynapses = numberOfSynapses + 1
                                                        end if
                                                    end do
                                                end if
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do
                end do
            end if  
            
            
            !Interneuron to MotorUnitPool
            if (size(interneuronPools)>0 .and. size(motorUnitPools)>0) then
                do poolOut = 1, size(interneuronPools)
                    do unitOut = 1, size(interneuronPools(poolOut)%unit)
                        neuralSource = trim(interneuronPools(poolOut)%pool) // '-' // &
                        interneuronPools(poolOut)%unit(unitOut)%neuronKind
                        interneuronPools(poolOut)%unit(unitOut)%SynapsesOut = conf%determineSynapses(neuralSource)
                        do synapseIn = 1, size(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item)
                            paramTag = 'Con:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,interneuronPools(poolOut)%pool, 0)
                            read(paramChar, *)conn
                            conn = conn / 100.0

                            paramTag = 'gmax:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,interneuronPools(poolOut)%pool, 0)
                            read(paramChar,*)gmax

                            paramTag = 'delay:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            
                            paramChar = conf%parameterSet(paramTag, interneuronPools(poolOut)%pool, 0)
                            read(paramChar, *)delay
                            
                            paramTag = 'dec:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, interneuronPools(poolOut)%pool, 0)
                            if (trim(paramChar)=='inf') then
                                declineFactor = 1e6
                            else 
                                read(paramChar, *)declineFactor
                            end if
                            
                            paramTag = 'dyn:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            dyn = conf%parameterSet(paramTag, interneuronPools(poolOut)%pool, 0)
                            
                            if (trim(dyn).ne.'None' .and. trim(dyn).ne.'none') then
                                paramTag = 'var:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                
                                paramChar = conf%parameterSet(paramTag,interneuronPools(poolOut)%pool, 0)
                                read(paramChar, *)var

                                paramTag = 'tau:' // trim(interneuronPools(poolOut)%pool) // '-' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(interneuronPools(poolOut)%&
                                            unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(interneuronPools(poolOut)%&
                                            unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%&
                                            item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,interneuronPools(poolOut)%pool, 0)
                                read(paramChar, *)tau
                            else
                                var = 0.0
                                tau = 1e6
                            end if
                            
                            do poolIn = 1, size(motorUnitPools)
                                if (trim(interneuronPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)==&
                                    trim(motorUnitPools(poolIn)%pool)) then
                                    do unitIn = 1, size(motorUnitPools(poolIn)%unit)
                                        do compartmentIn = 1, size(motorUnitPools(poolIn)%unit(unitIn)%Compartments)
                                            if ((trim(interneuronPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(1)%string)==&
                                                trim(motorUnitPools(poolIn)%pool)) .and.&
                                                (trim(interneuronPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(2)%string)==&
                                                trim(motorUnitPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                                (trim(interneuronPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(3)%string)==&
                                                trim(motorUnitPools(poolIn)%unit(unitIn)%&
                                                Compartments(compartmentIn)%compKind))) then
                                                call random_number(randomNumber)
                                                if (randomNumber <= conn) then
                                                    do synapseComp = 1, size(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                            Compartments(compartmentIn)%SynapsesIn)
                                                        if (motorUnitPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%synapseKind ==&
                                                            interneuronPools(poolOut)%unit(unitOut)%&
                                                            SynapsesOut%item(synapseIn)%item(4)%string) then
                                                            
                                                            if (declineFactor<1e5) then
                                                                neuronsDistance = abs(motorUnitPools(poolIn)%&
                                                                unit(unitIn)%position_mm &
                                                                            - interneuronPools(poolOut)%unit(unitOut)%position_mm)
                                                                weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                                gmax = gmax * weight
                                                            end if
                                                            call motorUnitPools(poolIn)%unit(unitIn)%&
                                                            Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                            addConductance(gmax, delay, dyn, var, tau)
                                                            
                                                            if (allocated(interneuronPools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)) then
                                                                
                                                                
                                                                allocate(tempTransmitSpikes(size(interneuronPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses)))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call tempTransmitSpikes(j)%&
                                                                    assignSynapse(interneuronPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                                end do
                                                                
                                                                
                                                                deallocate(interneuronPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses)
                                                                
                                                                allocate(interneuronPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call interneuronPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(j)%&
                                                                    assignSynapse(tempTransmitSpikes(j)%synapse)
                                                                end do

                                                                interneuronPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                        SynapsePointer()        

                                                                call interneuronPools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                                assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                Compartments(compartmentIn)%SynapsesIn(synapseComp))  
                                                                
                                                                deallocate(tempTransmitSpikes)
                                                            else 
                                                                allocate(interneuronPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1))
                                                                interneuronPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)=&
                                                                        SynapsePointer()        

                                                                call interneuronPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)%&
                                                                        assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                        Compartments(compartmentIn)%&
                                                                        SynapsesIn(synapseComp))
                                                            end if            
                                                            
                                                            newIndex = size(motorUnitPools(poolIn)%&
                                                                unit(unitIn)%Compartments(compartmentIn)%&
                                                                SynapsesIn(synapseComp)%gmax_muS)
                                                            
                                                            call integerAddToList(interneuronPools(poolOut)%&
                                                            unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                            
                                                            numberOfSynapses = numberOfSynapses + 1
                                                        end if
                                                    end do
                                                end if
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do
                end do
            end if       

            
            !MotorUnitPool to InterneuronPool
            if (size(motorUnitPools)> 0 .and. size(interneuronPools)>0) then
                do poolOut = 1, size(motorUnitPools)
                    do unitOut = 1, size(motorUnitPools(poolOut)%unit)
                        neuralSource = trim(motorUnitPools(poolOut)%pool) // '-' // motorUnitPools(poolOut)%unit(unitOut)%neuronKind
                        motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut = conf%determineSynapses(neuralSource)
                        do synapseIn = 1, size(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item)
                            paramTag = 'Con:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,motorUnitPools(poolOut)%pool, 0)
                            read(paramChar, *)conn
                            conn = conn / 100.0

                            paramTag = 'gmax:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,motorUnitPools(poolOut)%pool, 0)
                            read(paramChar,*)gmax

                            paramTag = 'delay:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            
                            paramChar = conf%parameterSet(paramTag, motorUnitPools(poolOut)%pool, 0)
                            read(paramChar, *)delay
                            
                            paramTag = 'dec:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, motorUnitPools(poolOut)%pool, 0)
                            if (trim(paramChar)=='inf') then
                                declineFactor = 1e6
                            else 
                                read(paramChar, *)declineFactor
                            end if
                            
                            paramTag = 'dyn:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            dyn = conf%parameterSet(paramTag, motorUnitPools(poolOut)%pool, 0)
                            
                            if (trim(dyn).ne.'None' .and. trim(dyn).ne.'none') then
                                paramTag = 'var:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                
                                paramChar = conf%parameterSet(paramTag,motorUnitPools(poolOut)%pool, 0)
                                read(paramChar, *)var

                                paramTag = 'tau:' // trim(motorUnitPools(poolOut)%pool) // '-' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(motorUnitPools(poolOut)%&
                                            unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,motorUnitPools(poolOut)%pool, 0)
                                read(paramChar, *)tau
                            else
                                var = 0.0
                                tau = 1e6
                            end if
                            
                            do poolIn = 1, size(interneuronPools)
                                if (trim(motorUnitPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)==&
                                    trim(interneuronPools(poolIn)%pool)) then
                                    do unitIn = 1, size(interneuronPools(poolIn)%unit)
                                        do compartmentIn = 1, size(interneuronPools(poolIn)%unit(unitIn)%Compartments)
                                            if ((trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(1)%string)==&
                                                trim(interneuronPools(poolIn)%pool)) .and.&
                                                (trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(2)%string)==&
                                                trim(interneuronPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                                (trim(motorUnitPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(3)%string)==&
                                                trim(interneuronPools(poolIn)%unit(unitIn)%&
                                                Compartments(compartmentIn)%compKind))) then
                                                call random_number(randomNumber)
                                                if (randomNumber <= conn) then
                                                    do synapseComp = 1, size(interneuronPools(poolIn)%unit(unitIn)%&
                                                                            Compartments(compartmentIn)%SynapsesIn)
                                                        if (interneuronPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%synapseKind ==&
                                                            motorUnitPools(poolOut)%unit(unitOut)%&
                                                            SynapsesOut%item(synapseIn)%item(4)%string) then
                                                            
                                                            if (declineFactor<1e5) then
                                                            neuronsDistance = abs(interneuronPools(poolIn)%unit(unitIn)%position_mm&
                                                                            - motorUnitPools(poolOut)%unit(unitOut)%position_mm)
                                                                weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                                gmax = gmax * weight
                                                            end if
                                                            call interneuronPools(poolIn)%unit(unitIn)%&
                                                            Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                            addConductance(gmax, delay, dyn, var, tau)
                                                            
                                                            if (allocated(motorUnitPools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)) then
                                                                
                                                                
                                                                allocate(tempTransmitSpikes(size(motorUnitPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses)))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call tempTransmitSpikes(j)%&
                                                                    assignSynapse(motorUnitPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                                end do                                                               
                                                                
                                                                deallocate(motorUnitPools(poolOut)%&
                                                                           unit(unitOut)%transmitSpikesThroughSynapses)
                                                                
                                                                allocate(motorUnitPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call motorUnitPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(j)%&
                                                                    assignSynapse(tempTransmitSpikes(j)%synapse)
                                                                end do

                                                                motorUnitPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                        SynapsePointer()        

                                                                call motorUnitPools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                                assignSynapse(interneuronPools(poolIn)%unit(unitIn)%&
                                                                Compartments(compartmentIn)%SynapsesIn(synapseComp))  
                                                                
                                                                deallocate(tempTransmitSpikes)
                                                            else 
                                                                allocate(motorUnitPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1))
                                                                motorUnitPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)=&
                                                                        SynapsePointer()        

                                                                call motorUnitPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)%&
                                                                        assignSynapse(interneuronPools(poolIn)%unit(unitIn)%&
                                                                        Compartments(compartmentIn)%&
                                                                        SynapsesIn(synapseComp))
                                                            end if            
                                                            
                                                            newIndex = size(interneuronPools(poolIn)%&
                                                                unit(unitIn)%Compartments(compartmentIn)%&
                                                                SynapsesIn(synapseComp)%gmax_muS)
                                                            
                                                            call integerAddToList(motorUnitPools(poolOut)%&
                                                            unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                            
                                                            numberOfSynapses = numberOfSynapses + 1
                                                        end if
                                                    end do
                                                end if
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do
                end do
            end if
            !Afferent to MotorUnitPool
            if (size(afferentPools)>0 .and. size(motorUnitPools)>0) then
                do poolOut = 1, size(afferentPools)
                    do unitOut = 1, size(afferentPools(poolOut)%unit)
                        neuralSource = trim(afferentPools(poolOut)%pool) // '-'&
                         // afferentPools(poolOut)%unit(unitOut)%neuronKind
                        afferentPools(poolOut)%unit(unitOut)%SynapsesOut = conf%determineSynapses(neuralSource)
                        do synapseIn = 1, size(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item)
                            paramTag = 'Con:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            
                            paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                            read(paramChar, *)conn
                            conn = conn / 100.0

                            paramTag = 'gmax:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)                            
                            paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                            read(paramChar,*)gmax

                            paramTag = 'delay:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)                            
                            paramChar = conf%parameterSet(paramTag, afferentPools(poolOut)%pool, 0)
                            read(paramChar, *)delay
                            
                            paramTag = 'dec:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, afferentPools(poolOut)%pool, 0)
                            if (trim(paramChar)=='inf') then
                                declineFactor = 1e6
                            else 
                                read(paramChar, *)declineFactor
                            end if
                            
                            paramTag = 'dyn:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            dyn = conf%parameterSet(paramTag, afferentPools(poolOut)%pool, 0)
                            
                            if (trim(dyn).ne.'None' .and. trim(dyn).ne.'none') then
                                paramTag = 'var:' // trim(afferentPools(poolOut)%pool) // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                
                                paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                                read(paramChar, *)var

                                paramTag = 'tau:' // trim(afferentPools(poolOut)%pool) // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                                SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                                read(paramChar, *)tau
                            else
                                var = 0.0
                                tau = 1e6
                            end if
                            
                            do poolIn = 1, size(motorUnitPools)
                                if (trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)==&
                                    trim(motorUnitPools(poolIn)%pool)) then
                                    do unitIn = 1, size(motorUnitPools(poolIn)%unit)
                                        do compartmentIn = 1, size(motorUnitPools(poolIn)%unit(unitIn)%Compartments)
                                            if ((trim(afferentPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(1)%string)==&
                                                trim(motorUnitPools(poolIn)%pool)) .and.&
                                                (trim(afferentPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(2)%string)==&
                                                trim(motorUnitPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                                (trim(afferentPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(3)%string)==&
                                                trim(motorUnitPools(poolIn)%unit(unitIn)%&
                                                Compartments(compartmentIn)%compKind))) then
                                                call random_number(randomNumber)
                                                if (randomNumber <= conn) then
                                                    do synapseComp = 1, size(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                            Compartments(compartmentIn)%SynapsesIn)
                                                        if (motorUnitPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%synapseKind ==&
                                                            afferentPools(poolOut)%unit(unitOut)%&
                                                            SynapsesOut%item(synapseIn)%item(4)%string) then
                                                            
                                                            if (declineFactor<1e5) then
                                                                neuronsDistance = abs(motorUnitPools(poolIn)%&
                                                                                unit(unitIn)%position_mm &
                                                                            - afferentPools(poolOut)%&
                                                                            unit(unitOut)%position_mm)
                                                                weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                                gmax = gmax * weight
                                                            end if
                                                            call motorUnitPools(poolIn)%unit(unitIn)%&
                                                            Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                            addConductance(gmax, delay, dyn, var, tau)
                                                            
                                                            if (allocated(afferentPools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)) then
                                                                
                                                                
                                                                allocate(tempTransmitSpikes(size(afferentPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses)))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call tempTransmitSpikes(j)%&
                                                                        assignSynapse(afferentPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                                end do
                                                                
                                                                
                                                                deallocate(afferentPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses)
                                                                
                                                                allocate(afferentPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call afferentPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(j)%&
                                                                    assignSynapse(tempTransmitSpikes(j)%synapse)
                                                                end do

                                                                afferentPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                        SynapsePointer()        

                                                                call afferentPools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                                assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                Compartments(compartmentIn)%SynapsesIn(synapseComp))  
                                                                
                                                                deallocate(tempTransmitSpikes)
                                                            else 
                                                                allocate(afferentPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1))
                                                                afferentPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)=&
                                                                        SynapsePointer()        

                                                                call afferentPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)%&
                                                                        assignSynapse(motorUnitPools(poolIn)%unit(unitIn)%&
                                                                        Compartments(compartmentIn)%&
                                                                        SynapsesIn(synapseComp))
                                                            end if            
                                                            
                                                            newIndex = size(motorUnitPools(poolIn)%&
                                                                unit(unitIn)%Compartments(compartmentIn)%&
                                                                SynapsesIn(synapseComp)%gmax_muS)
                                                            
                                                            call integerAddToList(afferentPools(poolOut)%&
                                                            unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                            
                                                            numberOfSynapses = numberOfSynapses + 1
                                                        end if
                                                    end do
                                                end if
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do
                end do
            end if  

            !Afferent to InterNeuronPool
            if (size(afferentPools)>0 .and. size(interneuronPools)>0) then
                do poolOut = 1, size(afferentPools)
                    do unitOut = 1, size(afferentPools(poolOut)%unit)
                        neuralSource = trim(afferentPools(poolOut)%pool) // '-'&
                         // afferentPools(poolOut)%unit(unitOut)%neuronKind
                        afferentPools(poolOut)%unit(unitOut)%SynapsesOut = conf%determineSynapses(neuralSource)
                        do synapseIn = 1, size(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item)
                            paramTag = 'Con:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                            read(paramChar, *)conn
                            conn = conn / 100.0

                            paramTag = 'gmax:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)                            
                            paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                            read(paramChar,*)gmax

                            paramTag = 'delay:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)                            
                            paramChar = conf%parameterSet(paramTag, afferentPools(poolOut)%pool, 0)
                            read(paramChar, *)delay
                            
                            paramTag = 'dec:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, afferentPools(poolOut)%pool, 0)
                            if (trim(paramChar)=='inf') then
                                declineFactor = 1e6
                            else 
                                read(paramChar, *)declineFactor
                            end if
                            
                            paramTag = 'dyn:' // trim(afferentPools(poolOut)%pool) // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(4)%string)
                            dyn = conf%parameterSet(paramTag, afferentPools(poolOut)%pool, 0)
                            
                            if (trim(dyn).ne.'None' .and. trim(dyn).ne.'none') then
                                paramTag = 'var:' // trim(afferentPools(poolOut)%pool) // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                
                                paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                                read(paramChar, *)var

                                paramTag = 'tau:' // trim(afferentPools(poolOut)%pool) // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%neuronKind) // '>'&
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                                SynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(afferentPools(poolOut)%unit(unitOut)%&
                                            SynapsesOut%item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,afferentPools(poolOut)%pool, 0)
                                read(paramChar, *)tau
                            else
                                var = 0.0
                                tau = 1e6
                            end if
                            
                            do poolIn = 1, size(interneuronPools)
                                if (trim(afferentPools(poolOut)%unit(unitOut)%SynapsesOut%item(synapseIn)%item(1)%string)==&
                                    trim(interneuronPools(poolIn)%pool)) then
                                    do unitIn = 1, size(interneuronPools(poolIn)%unit)
                                        do compartmentIn = 1, size(interneuronPools(poolIn)%unit(unitIn)%Compartments)
                                            if ((trim(afferentPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(1)%string)==&
                                                trim(interneuronPools(poolIn)%pool)) .and.&
                                                (trim(afferentPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(2)%string)==&
                                                trim(interneuronPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                                (trim(afferentPools(poolOut)%unit(unitOut)%&
                                                    SynapsesOut%item(synapseIn)%item(3)%string)==&
                                                trim(interneuronPools(poolIn)%unit(unitIn)%&
                                                Compartments(compartmentIn)%compKind))) then
                                                call random_number(randomNumber)
                                                if (randomNumber <= conn) then
                                                    do synapseComp = 1, size(interneuronPools(poolIn)%unit(unitIn)%&
                                                                            Compartments(compartmentIn)%SynapsesIn)
                                                        if (interneuronPools(poolIn)%unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%synapseKind ==&
                                                            afferentPools(poolOut)%unit(unitOut)%&
                                                            SynapsesOut%item(synapseIn)%item(4)%string) then
                                                            
                                                            if (declineFactor<1e5) then
                                                                neuronsDistance = abs(interneuronPools(poolIn)%&
                                                                                unit(unitIn)%position_mm &
                                                                                - afferentPools(poolOut)%&
                                                                                unit(unitOut)%position_mm)
                                                                weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                                gmax = gmax * weight
                                                            end if
                                                            call interneuronPools(poolIn)%unit(unitIn)%&
                                                            Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                            addConductance(gmax, delay, dyn, var, tau)
                                                            
                                                            if (allocated(afferentPools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)) then                                                                
                                                                
                                                                allocate(tempTransmitSpikes(size(afferentPools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses)))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call tempTransmitSpikes(j)%&
                                                                        assignSynapse(afferentPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                                end do
                                                                
                                                                
                                                                deallocate(afferentPools(poolOut)%&
                                                                        unit(unitOut)%transmitSpikesThroughSynapses)
                                                                
                                                                allocate(afferentPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                                
                                                                do j = 1, size(tempTransmitSpikes)
                                                                    call afferentPools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(j)%&
                                                                    assignSynapse(tempTransmitSpikes(j)%synapse)
                                                                end do

                                                                afferentPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                        SynapsePointer()        

                                                                call afferentPools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                                assignSynapse(interneuronPools(poolIn)%unit(unitIn)%&
                                                                Compartments(compartmentIn)%SynapsesIn(synapseComp))  
                                                                
                                                                deallocate(tempTransmitSpikes)
                                                            else 
                                                                allocate(afferentPools(poolOut)%&
                                                                        unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1))
                                                                afferentPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)=&
                                                                        SynapsePointer()        

                                                                call afferentPools(poolOut)%unit(unitOut)%&
                                                                        transmitSpikesThroughSynapses(1)%&
                                                                        assignSynapse(interneuronPools(poolIn)%unit(unitIn)%&
                                                                        Compartments(compartmentIn)%&
                                                                        SynapsesIn(synapseComp))
                                                            end if            
                                                            
                                                            newIndex = size(interneuronPools(poolIn)%&
                                                                unit(unitIn)%Compartments(compartmentIn)%&
                                                                SynapsesIn(synapseComp)%gmax_muS)
                                                            
                                                            call integerAddToList(afferentPools(poolOut)%&
                                                            unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                            
                                                            numberOfSynapses = numberOfSynapses + 1
                                                        end if
                                                    end do
                                                end if
                                            end if
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do
                end do
            end if  

            print '(A, I10, A)', 'All the ', numberOfSynapses,  ' synapses were built'

            numberOfSynapticNoise = 0
            neuralSource = 'Noise'            
            NoiseSynapsesOut = conf%determineSynapses(neuralSource)    

            if (allocated(synapticNoisePools)) deallocate(synapticNoisePools)
            if (allocated(NoiseSynapsesOut%item) .and. size(interneuronPools)>0) then
                do synapseIn = 1, size(NoiseSynapsesOut%item)
                    do poolIn = 1, size(interneuronPools)
                        if (trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string)==&
                            trim(interneuronPools(poolIn)%pool)) then
                            if (allocated(synapticNoisePools)) then
                                allocate(tempSynNoise(size(synapticNoisePools)))
                                do i = 1, size(synapticNoisePools)
                                    tempSynNoise(i) = synapticNoisePools(i)
                                end do
                                deallocate(synapticNoisePools)
                                allocate(synapticNoisePools(size(tempSynNoise)+1))
                                do i = 1, size(tempSynNoise)
                                    synapticNoisePools(i) = tempSynNoise(i)
                                end do
                                synapticNoisePools(size(tempSynNoise)+1) = &
                                SynapticNoise(conf, NoiseSynapsesOut%item(synapseIn)%item(1)%string)
                                deallocate(tempSynNoise)
                            else 
                                allocate(synapticNoisePools(1))
                                synapticNoisePools(1) = SynapticNoise(conf, NoiseSynapsesOut%item(synapseIn)%item(1)%string)
                            end if
                            poolOut = size(synapticNoisePools)
                            paramTag = 'gmax:Noise>' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag,synapticNoisePools(poolOut)%pool, 0)
                            read(paramChar,*)gmax

                            paramTag = 'delay:Noise>' &
                                    // trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string) &
                                    // '-' &
                                    // trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string) &
                                    // '@' &
                                    // trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string) &
                                    // '|' &
                                    // trim(NoiseSynapsesOut%item(synapseIn)%item(4)%string)
                                    
                            paramChar = conf%parameterSet(paramTag, synapticNoisePools(poolOut)%pool, 0)
                            read(paramChar, *)delay
                                    
                            paramTag = 'dec:Noise>' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(4)%string)
                            paramChar = conf%parameterSet(paramTag, synapticNoisePools(poolOut)%pool, 0)
                            if (trim(paramChar)=='inf') then
                                declineFactor = 1e6
                            else 
                                read(paramChar, *)declineFactor
                            end if
                                
                            paramTag = 'dyn:Noise>'&
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string) &
                                        // '-' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string) &
                                        // '@' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string) &
                                        // '|' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(4)%string)
                            dyn = conf%parameterSet(paramTag, synapticNoisePools(poolOut)%pool, 0)
                                    
                            if (trim(dyn).ne.'None' .and. trim(dyn).ne.'none') then
                                paramTag = 'var:Noise>'&
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string)&
                                        // '-' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string)&
                                        // '@' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string)&
                                        // '|' &
                                        // trim(NoiseSynapsesOut%item(synapseIn)%item(4)%string)
                                        
                                paramChar = conf%parameterSet(paramTag,synapticNoisePools(poolOut)%pool, 0)
                                read(paramChar, *)var

                                paramTag = 'tau:Noise>'&
                                            // trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string)&
                                            // '-' &
                                            // trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string)&
                                            // '@' &
                                            // trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string)&
                                            // '|' &
                                            // trim(NoiseSynapsesOut%item(synapseIn)%item(4)%string)
                                paramChar = conf%parameterSet(paramTag,synapticNoisePools(poolOut)%pool, 0)
                                read(paramChar, *)tau
                            else
                                var = 0.0
                                tau = 1e6
                            end if
                            
                            do unitOut = 1, size(synapticNoisePools(poolOut)%unit)
                                    unitIn = unitOut
                                    do compartmentIn = 1, size(interneuronPools(poolIn)%unit(unitIn)%Compartments)
                                        if ((trim(NoiseSynapsesOut%item(synapseIn)%item(1)%string)==&
                                            trim(interneuronPools(poolIn)%pool)) .and.&
                                            (trim(NoiseSynapsesOut%item(synapseIn)%item(2)%string)==&
                                            trim(interneuronPools(poolIn)%unit(unitIn)%neuronKind)) .and.&
                                            (trim(NoiseSynapsesOut%item(synapseIn)%item(3)%string)==&
                                            trim(interneuronPools(poolIn)%unit(unitIn)%&
                                                    Compartments(compartmentIn)%compKind))) then
                                            
                                            do synapseComp = 1, size(interneuronPools(poolIn)%unit(unitIn)%&
                                                                Compartments(compartmentIn)%SynapsesIn)                                                    
                                                    if (NoiseSynapsesOut%item(synapseIn)%item(4)%string==&
                                                        interneuronPools(poolIn)%unit(unitIn)%&
                                                        Compartments(compartmentIn)%SynapsesIn(synapseComp)%synapseKind) then
                                                        if (declineFactor<1e5) then
                                                            neuronsDistance = abs(interneuronPools(poolIn)%&
                                                                            unit(unitIn)%position_mm&
                                                                            - synapticNoisePools(poolOut)%unit(unitOut)%&
                                                                            position_mm)
                                                            weight = declineFactor / (declineFactor + neuronsDistance**2)
                                                            gmax = gmax * weight
                                                        end if
                                                        call interneuronPools(poolIn)%unit(unitIn)%&
                                                        Compartments(compartmentIn)%SynapsesIn(synapseComp)%&
                                                        addConductance(gmax, delay, dyn, var, tau)
                                                        
                                                        if (allocated(synapticNoisePools(poolOut)%&
                                                            unit(unitOut)%transmitSpikesThroughSynapses)) then
                                                            
                                                            
                                                            allocate(tempTransmitSpikes(size(synapticNoisePools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses)))
                                                            
                                                            do j = 1, size(tempTransmitSpikes)
                                                                call tempTransmitSpikes(j)%&
                                                                assignSynapse(synapticNoisePools(poolOut)%&
                                                                unit(unitOut)%transmitSpikesThroughSynapses(j)%synapse)
                                                            end do
                                                            
                                                            
                                                            deallocate(synapticNoisePools(poolOut)%&
                                                                    unit(unitOut)%transmitSpikesThroughSynapses)
                                                            
                                                            allocate(synapticNoisePools(poolOut)%&
                                                                    unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1))
                                                            
                                                            do j = 1, size(tempTransmitSpikes)
                                                                call synapticNoisePools(poolOut)%unit(unitOut)%&
                                                                transmitSpikesThroughSynapses(j)%&
                                                                assignSynapse(tempTransmitSpikes(j)%synapse)
                                                            end do

                                                            synapticNoisePools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)=&
                                                                    SynapsePointer()        

                                                            call synapticNoisePools(poolOut)%unit(unitOut)%&
                                                            transmitSpikesThroughSynapses(size(tempTransmitSpikes)+1)%&
                                                            assignSynapse(interneuronPools(poolIn)%unit(unitIn)%&
                                                            Compartments(compartmentIn)%SynapsesIn(synapseComp))
                                                            
                                                            deallocate(tempTransmitSpikes)
                                                        else
                                                            allocate(synapticNoisePools(poolOut)%&
                                                                    unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(1))
                                                            synapticNoisePools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(1)=&
                                                                    SynapsePointer()

                                                            call synapticNoisePools(poolOut)%unit(unitOut)%&
                                                                    transmitSpikesThroughSynapses(1)%&
                                                                    assignSynapse(interneuronPools(poolIn)%unit(unitIn)%&
                                                                    Compartments(compartmentIn)%&
                                                                    SynapsesIn(synapseComp))
                                                        end if            
                                                        
                                                        newIndex = size(interneuronPools(poolIn)%&
                                                            unit(unitIn)%Compartments(compartmentIn)%&
                                                            SynapsesIn(synapseComp)%gmax_muS)
                                                        
                                                        call integerAddToList(synapticNoisePools(poolOut)%&
                                                        unit(unitOut)%indicesOfSynapsesOnTarget,newIndex)
                                                        numberOfSynapticNoise = numberOfSynapticNoise + 1
                                                    end if
                                                end do                                            
                                            end if
                                    end do                                    
                            end do
                        end if
                    end do
                end do
            end if

            if (.not. allocated(synapticNoisePools)) allocate(synapticNoisePools(0))

            if (allocated(tempSynNoise)) deallocate(tempSynNoise)
            if (allocated(tempTransmitSpikes)) deallocate(tempTransmitSpikes)

            print '(A, I10, A)', 'All the ', numberOfSynapticNoise,  ' synaptic noises were built'


        end function

end module SynapsesFactoryModule 