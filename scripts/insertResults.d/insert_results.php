<?php

function insert_results( $STAARfilename )
{
  /* Connect to database */

  /* Read the results file into RAM */
  $STAARfile = file($STAARfilename, FILE_IGNORE_NEW_LINES);

  $previous_PDB = "";

  /* Parse each line of the file */
  foreach ($STAARfile as $line_num=>$line)
    {
      /* Skip comment lines */
      if($line[0] == '#') continue;
      $fields = explode(',', $line);

      /* Check to see if the PDB exists in the database already */
      $PDB = explode('/', $fields[9]);
      $PDB = explode('.', $PDB[count($PDB)-1]);
      if($PDB != $prevPDB)
        {
          $prevPDB = $PDB;
          //print "SELECT id FROM pdb where pdb=\"".$PDB[0]."\"\n";
          //$results = mysql_query("SELECT id FROM pdb where pdb=\"".$PDB[0]."\"");
      
          /* If it doesn't exist already, we need in insert the PDB id 
             into the database first                                   */
          if( mysql_num_rows($result) != 1 )
            {
              $query = "INSERT INTO pdb (pdb) VALUES(\"".$PDB[0]."\")";
              $result = mysql_query($query);
              if( !$result ) die('Invalid query: '.mysql_error());
              $PDBid = mysql_insert_id();
            }
          else
            {
              $PDBid = mysql_fetch_array($result, MYSQL_NUM);
              $PDBid = $PDBid[0];
            }
        }

      /* Now, the definitely PDB exists, let's insert the results */
      
      
    }

}

insert_results("/lustre/AQ/Ghost/STAARwEnergy.csv");

?>
