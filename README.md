# MSRB-Fagmenter's README
## Compiling MSRB-Fragmenter with Eclipse
<p>This tutorial assumes that you have Eclipse installed, as well as the Maven Development Tools (m2e 1.6+) 0.30. If you haven't please select the section "Setting up Eclipse and Maven". The current version of the msrb-fragemnter requires Java 8.</p>

1) Cloning the project repository directly from Eclipse
> -  Click on import projects
> - Select Git
> - Fill out the form: Usually, adding the https URI will prompt an automated filling of the other fields (repository path, host name, username)
> - add the user password
> - Click on ''
> - In the Destination menu, point to the correct destination
> - Select the correct branch you want
> - Then select 'Next'
> - In the window 'Select a wizard to use for importing projects', selectt 'Import existing Eclipse projects' and click 'Next'
> - In the "Import Projects" window, select the correct project, and click 'Finish'


2) Import the repository From an existing local project , if is has been cloned locally
> - Select Import projects -> General -> Existing Projects into Workspace -> Select root directory

3) If you get an error message for missing src folders, e.g. "The project cfmid is missing source folder '/src/main/resources' "
> - Right-Click on the Project Name
> - Select Build Path -> Configure Build Path
> - In the appearing window, click on the 'source' button:
<ul>
<li>Select each of the missing folders, and clock on the 'Remove' button on the right.</li>
<li>Select 'Apply'</li>
<li>On the right menu: click "Add Folder" -> "Create Folder"</li>
<li>Add the path of the missing folder (e.g. /src/main/resources) under 'Folder name:'</li>
<li>Click 'Finish'</li>
<li>Repeat this for any missing folder</li>
</ul>

4) Selecting the correct Java version for the project
> - Go to the Build Path -> Configure Build Path
> - Click on 'Libraries', then select the current version of the JRE to be removed, and click on the "Remove" button.
> - Then select "Add Library" -> "JRE System Library"
> - On the new window, select the correct JRE version

5) Creating an executable
> - Right-click on the project's name and select "export"
> - In the newly opened window, select Java -> Runnable JAR file
> - Select Next
> - In the "Runnable JAR File export" window:
	<ul>
	<li>Select "RuleBasedFrag" for the launch configuration</li>
	<li>Select "cfmid_plus/bin/msrb-fragmenter.jar" as Export destination</li>
	<li>Make sure to select "Extract required libraries into generated JAR" in the Library handling section</li>
	<li>Check the box "Save as ANT script", and add the ANT script location</li>
	<li>Click Finish</li>
	</ul>