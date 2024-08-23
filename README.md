# Reconstructor
This tool is for the curation of reconstruction reactions. It allows users to preview, save and run tests on reactions of their choice.

It is live on http://reconstructor.chatimd.org/

## Installation and Setup (to run locally)
To get started with the Reconstruction Curation Tool, follow these steps to set up the environment on your local machine.
### Prerequisites
Ensure you have the following installed:
- Python (version 3.10.12)
- Django (version 5.0.1)
- RDT - https://github.com/asad/ReactionDecoder/releases (v2.4.1)

### Run the Application
1. Make sure you have postgresql or another DB management system
2. Specify the database name, user and password in reactions_project/settings.py\
3. Execute the following commands:\
`pip install -r requirements. txt`\
`python manage.py runserver`
## Code Structure
### HTML Files
1. **Home_page.html**
    - Main page for inputting reaction data.
    - Includes form fields for substrates, products, and reaction direction.
    - Integrates Django template language for dynamic content rendering.
    - Contains modals for saving reactions, adding information, and ChemDoodle sketcher.
    - Links to CSS for styling and JavaScript for interactive elements.
2. **saved_reactions.html**
    - Displays a list of saved reactions for a given user.
    - Each reaction entry includes details like name, substrates, products, and chemical formula.
    - Provides options to view more details or delete a reaction. **(this is now in static/reactions/js_savedreactions/reactionViewAndDelete.js)**
    - Add to VMH Button Handling: Integrates a user interface component for submitting selected reactions to be added to the VMH database **(this is now in static/reactions/js_savedreactions/hadleAdd2VMH.js)** :
        - JavaScript functionality to manage user selections of reactions for submission.
        - Dynamic generation of a modal window that previews selected reactions, allowing users to review and provide additional details before submission.
        - Validation to ensure all required information is provided for each selected reaction, with user-friendly error messages for missing data.
        - Communication with the backend views.py function add_to_vmh through AJAX, submitting the detailed reaction information.
        - Receive feedback on the process, including loading indicators and success or error messages upon completion.
3. **chemdoodle_sketcher.html**
    - Dedicated page for the ChemDoodle sketcher.
    - Includes necessary ChemDoodle Web components and scripts.
    - Provides functions for getting molecule data, loading molecules into the sketcher, and clearing the sketcher.
    - This is used as a sub-page in the **Home_page.html** when a user wants to draw a molecule.
### CSS Files
1. **input_reaction.css**
    - Styles for the **Home_page.html** page.
    - Defines layout for form fields, buttons, modals, and containers.
    - Ensures responsive design for different screen sizes.
    - Includes custom styles for error messages and status indicators.
2. **saved_reactions.css**
    - Styles for the **saved_reactions**.html page.
    - Focuses on list styling and reaction information presentation.
    - Customizes buttons for viewing and deleting reactions.
    - Ensures readability and clear layout for reaction lists.
### JavaScript Files
1. **chemdoodleIntegration.js**
    - Modular Integration with ChemDoodle: Provides functions to integrate and manage ChemDoodle sketcher within the application.
    - Modal Display Control: Contains `hideChemdoodlestatus` function for showing or hiding the ChemDoodle sketcher modal based on user interaction.
    - Event Listener Attachment: Implements `attachEventListenersToSelects` for attaching event listeners to select elements, allowing dynamic interaction with the sketcher.
    - Dynamic Button Management: Functions `addStartDrawingButton`, `removeStartDrawingButton`, and `removeEditDrawingButton` manage the addition and removal of drawing-related buttons based on user selections.
    - ChemDoodle Interaction: `setCurrentlyDrawing` and `clearCurrentlyDrawing` functions manage the state of which molecule is currently being drawn or edited.
    - Molecule Data Handling: `saveDrawing` function handles the extraction of molecule data from the ChemDoodle sketcher, prompts for a molecule name, and updates the corresponding input fields.
    - Editing Functionality: `editDrawing` and `replaceStartWithEditButton` provide functionality to edit existing molecule drawings, including retrieving and sending data back to the ChemDoodle sketcher.
2. **formHandling.js**
    - Dynamic Form Submission and Handling: Manages the submission of various forms within the application, including the addition of reaction information, reaction form submission, and saving reactions.
    - Event Listener for Adding Information: Implements an event listener on the 'Add Information' button. Captures and sends user input data such as user key, information type, text, and external link type using AJAX POST request to `addInfo2Reaction` endpoint.
    - Error Handling and Response Management: Processes server responses, handling both success and error cases, and updates the user interface accordingly.
    - Reaction Form Submission Logic: Attaches an event listener to the main reaction form to handle submission. Prevents default form submission, manages form data, and sends it to the inputReactionUrl endpoint.
    - Loading Indicator Management: Displays a loading indicator during data processing and disables the submit button to prevent multiple submissions.
    - Error Display for Reaction Submission: Shows an error message in the user interface if an unexpected error occurs during form submission.
    - VMH (Virtual Metabolic Human) Integration: Manages the functionality to retrieve reaction data from VMH. It includes enabling the 'Get from VMH' button based on input and handling the button click event to fetch data from the `getFromVmhUrl` endpoint.
    - View Saved Reactions Feature: Attaches an event listener to the 'View Saved Reactions' button, prompting the user to enter a key and then navigating to the saved reactions page with the provided key.
    - Save Reaction Functionality: Handles the saving of a reaction. This involves capturing user input for reaction identification and user key, sending the data to the `saveReaction` endpoint, and providing feedback to the user on success or failure.
    - UX Enhancements: Includes various user experience improvements such as clearing input fields upon successful data submission, handling form field validation, and providing immediate feedback to user actions.
3. **inputFieldHandling.js**
    - Dynamic Field Management: Specializes in dynamically adding and removing input fields for substrates and products in the form.
    - Adding Substrate and Product Fields: Implements `addField` function, called when 'Add Substrate' or 'Add Product' buttons are clicked. This function creates a new set of input fields including a number input for stoichiometry, a text input for the compound, and a select input for the compound type.
    - Remove Button Functionality: Integrates event listeners on dynamically created 'Remove' buttons to delete their respective input group.
    - File Input Management: Handles the toggling of file input fields with toggleFileInput function, which switches between text and file inputs based on the selected compound type, specifically for 'MDL Mol file' type.
    - Field Creation with Pre-filled Data: Contains `addFieldWithData` function for adding fields with existing data, useful for pre-populating the form when editing existing reactions or loading data from external sources.
    - Updating Form Fields from External Data: Features `updateFormFields` function to clear existing input groups and create new ones based on provided data, particularly useful for loading reaction details from an external source like VMH (Virtual Metabolic Human).
    - Enhanced User Interaction: Improves user experience by allowing users to dynamically manage the number of substrates and products in their reactions and by providing a clear and intuitive interface for inputting or editing reaction components.
4. **inputFieldHandling.js**
    - Modal Display Control: Manages the display and hiding of modals in the application, specifically for 'Add Information' and 'Save Reaction' functionalities.
    'Add Information' Modal Handling:
        - Implements functionality to show the 'Add Information' modal upon clicking the 'Add Info' button.
        - Clears text content and response messages in the modal when opened.
        - Closes the modal when the close button is clicked.
    - Dynamic Placeholder Management:
        - Updates the placeholder text in the information text input field based on the selected information type (Reference, External Link, Gene Info, Comment).
        - Dynamically adds or removes an additional dropdown for the 'External Link' type, enabling further specificity for the link type (e.g., KEGG orthology, Wikipedia).
    - 'Save Reaction' Modal Handling:
        - Handles the display of the 'Save Reaction' modal when the 'Save Reaction' button is clicked.
        - Includes functionality to close the modal when the close button is clicked or when a click is detected outside the modal content.
    - Enhanced User Interaction:
        - Provides a user-friendly interface for adding reaction information and saving reactions, with clear and intuitive modal interactions.
        - Ensures a smooth and non-disruptive user experience by handling modal interactions efficiently and effectively.
5. **reactionDisplay.js**
    - Tab Management for Reaction Details: Handles the display logic for different tabs (Chem Info, References, External Links, Gene Info, Comments) in the reaction detail section.
        - Utilizes radio buttons to toggle between different content tabs.
        - Adds event listeners to each tab to show the corresponding content when selected.
    - Display Reaction Details: Implements `displayReactionDetails` function to show reaction data fetched from the server.
        - Handles display of error messages and links to external resources (like VMH URLs) when available.
        - Updates UI elements like status dots based on the data received (e.g., indicating whether substrates and products were found in databases).
        - Fetches and displays additional details for the reaction in various tabs.
    - Reaction Data Presentation: Contains `displayReactionData` function for visual representation of reaction details.
        - Dynamically generates sections for each reaction, including images, molecular formulas, and balancing information.
        - Creates and populates div elements with atom lists and other reaction-related data.
    - Atom List Population: Features `populateAtomList` function to display lists of atoms involved in substrates or products, including their counts and names.
        - Sorts atoms based on their count and displays them in a structured format.
    - Tab Content Display: Provides `displayTabContent` function to populate specific tabs (like References or External Links) with relevant data in a table format.
        -Generates HTML table elements dynamically based on the type of data received for each tab.
    - Status Dot Updates: Implements `updateStatusDots` function to update the visual indicators (dots) next to substrates and products, showing whether they are found in the database.
        - Makes these dots clickable, linking to external resources if available.
    - Response Message Updates: Includes `updateResponseMessage` function to display success or error messages based on the responses received from server interactions.
        - Dynamically updates the class of the message container to reflect the nature of the message (success or error).
6. **reactionInitialization.js**
    - Initial Reaction Data Loading: For cases where a specific reaction ID is passed in the URL, this initializes the reaction details when the page is loaded.
    - URL Parameter Handling:
        - Extracts the reaction_id from the URL query parameters using `URLSearchParams`.
        - Ensures that the application responds dynamically to URL changes, allowing for deep linking to specific reactions.
    - Fetching Reaction Data:
        - If a `reaction_id` is present in the URL, makes a fetch request to the server using the `getReaction` URL concatenated with the `reaction_id`.
        - Processes the response to display the reaction details on the page.
### Django Application Files
1. **views.py**
    - Django View Functions: Contains various view functions to handle different aspects of the web application.
    - Chemdoodle Sketcher: Renders the ChemDoodle sketcher page.
    - Input Reaction: Processes the POST request for a reaction input form, validating and saving the reaction data.
    - Fetch Reaction Details: Retrieves specific reaction details by ID and returns them in JSON format.
    - Add Information to Reaction: Adds additional information (like references, external links, gene info, comments) to a specific reaction.
    - Get Reaction Details: Fetches detailed information for a specific reaction, including references and external links.
    - User Key Validation: Validates a user key for authenticity.
    - Save User Reaction: Saves a user's reaction with additional details like a short name.
    - Saved Reactions: Renders a page listing all reactions saved by a specific user.
    - Delete Reaction: Handles the deletion of a saved reaction from a user's list.
    - Add to VMH: Handles the addition of new reactions and their associated metabolites to the Virtual Metabolic Human (VMH) database. This process involves:
        - Validating the user's permission to add data to the VMH.
        - Extracting and processing reaction details from the request.
        - Generating unique abbreviations for new metabolites and reactions not currently in the VMH database.
        - Saving updated reaction information with user-provided references, external links, and comments.
        - Adding new metabolites to the database (if needed), including their formulas, charges, and InChIKeys, using MATLAB integration.
        - Constructing and saving new reaction entities with complete details including subsystem, directionality, and associated gene information.
        - Providing feedback on the operation's success, including detailed lists of added reactions and metabolites.
2. **reaction_info.py**
    - Reaction Information Processing: Provides functions for processing chemical reactions, like balancing checks and molecular formula generation.
    - Total Charge Calculation: Calculates the total charge of molecules in a reaction.
    - Element Counting: Counts the elements in a molecule.
    - Reaction Balancing Checks: Determines whether a reaction is balanced in terms of atom count and charge.
    - Molecular Formula Generation: Constructs molecular formulas for reactions.
3. **models.py**
    - User Model: Represents users with unique keys and optional names. Includes a many-to-many relationship with the Reaction model to store saved reactions.
    - Reaction Model: Stores detailed information about chemical reactions.
        - Includes fields for substrates, products, their atom counts, charges, and molecular formulas.
        - Contains fields for storing JSON strings, such as subs_found, subs_miriams, prod_found, and prod_miriams, which hold information about whether substrates/products were found in databases and their corresponding identifiers.
        - VMH (Virtual Metabolic Human) specific fields (vmh_found, vmh_url, vmh_formula) indicate whether the reaction is found in the VMH database and store related data.
        - Information fields (references, ext_links, gene_info, comments) use JSONField to store various types of additional information related to the reaction.
4. **admin.py**
    - UserAdmin Class: Configures the admin interface for the User model.
        - list_display: Specifies the fields of the User model to display in the Django admin list view.
        - search_fields: Defines the fields on which admin users can search.
        - filter_horizontal: Enhances the user interface for editing many-to-many relationships (saved reactions for users).
    - ReactionAdmin Class: Configures the admin interface for the Reaction model.
        - list_display: Specifies the fields of the Reaction model to display in the Django admin list view.
        - search_fields: Defines the fields on which admin users can search.
    - Admin Site Registration: Registers the User and Reaction models along with their respective admin classes to the Django admin site, enabling admins to manage these models through the admin interface.
5. **utils/RDT.py**
    - Reaction Decoder Tool (RDT) Processing: Handles the use of RDT for generating atom-atom mapping and visualizing chemical reactions.
    - Placeholder Replacement: Replaces placeholders in the reaction file with actual values for labeling.
6. **utils/to_smiles.py**
    - Conversion to SMILES Strings: Converts various types of molecular identifiers to SMILES strings.
        - Supports conversions from VMH, SwissLipids, MetaNetX, CHeBI, and MDL Mol files.
    - Explicit Hydrogens in SMILES: Ensures all hydrogen atoms are explicitly represented in SMILES strings.
7. **utils/search_vmh.py**
    - VMH Database Interactions: Includes functions for handling chemical compounds and reactions with the VMH database.
    - Molecule Search: Searches for molecules in VMH and fetches MIRIAM IDs.
    - Reaction Check in VMH: Checks if a chemical reaction exists in the VMH database.
8. **utils/add_to_vmh_utils.py**
    - Gather Reaction Details: Extracts additional information for reactions such as direction, subsystems, and associated gene information.
    - JSON Preparation for MATLAB: Creates temporary JSON files containing metabolite and reaction details for MATLAB processing.
    - MATLAB Integration: Executes MATLAB scripts to add metabolites and reactions to the VMH, leveraging the COBRA Toolbox.
    - SMILES Conversion: Includes utilities for converting molecular identifiers into SMILES strings, generating InChIKeys, and calculating formulas and charges.
    - Metabolite and Reaction Abbreviation Generation: Generates unique abbreviations for new metabolites and reactions, facilitating their addition to the VMH database.
    - Error Handling and Feedback: Provides detailed error messages and success feedback, including the generation of unique identifiers for successfully added entities.