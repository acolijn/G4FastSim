#include "EventAction.hh"
#include "EventActionMessenger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

/**
 * @class EventActionMessenger
 * 
 * @brief Messenger class for EventAction.
 * 
 * This class is responsible for handling user commands related to EventAction.
 * It provides commands to set the spatial and time thresholds for clustering.
 */
namespace G4Sim {

/**
 * @brief Constructs an EventActionMessenger object.
 * 
 * @param eventAction Pointer to the EventAction object.
 */
EventActionMessenger::EventActionMessenger(EventAction* eventAction)
    : G4UImessenger(), fEventAction(eventAction) {
    
    // Command to set the spatial threshold


EventActionMessenger::~EventActionMessenger() {
}

/**
 * @brief Sets the new value for a given command.
 * 
 * This method is called when a new value is set for a command in the EventActionMessenger class.
 * It checks the command and updates the corresponding value in the EventAction class.
 * 
 * @param command The G4UIcommand object representing the command.
 * @param newValue The new value to be set.
 */
void EventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

}

} // namespace G4XamsSim
