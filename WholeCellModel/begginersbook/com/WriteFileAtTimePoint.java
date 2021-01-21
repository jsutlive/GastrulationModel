package begginersbook.com;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class WriteFileAtTimePoint
{
    public static void main(String[] args)
    {
        BufferedWriter bw = null;
        WriteToFile(bw, null, null);
    }
// Change file path to fit your computer's file system
    public static void WriteToFile(BufferedWriter bw, String myFileName, String myContent) {
        String path = "";
        if (myContent == null)
        {
            myContent = "This string would be written" +
                    " to the specified File";
        }
        if(myFileName == null)
        {
            path = "I:/Documents/Dartmouth/Code/myFile.txt";
        }
        else
        {
            path = "I:/Documents/Dartmouth/Code/DF10/" + myFileName + ".txt";
        }
        try
        {
            System.out.println("PATH: " + path);
            //specify the file name and path here:
            File file = new File(path);

            if(!file.exists()){file.createNewFile();}
            FileWriter fw = new FileWriter(file);
            bw = new BufferedWriter(fw);
            bw.write(myContent);
            System.out.println("File written Successfully");
        } catch (IOException ioe)
        {
            ioe.printStackTrace();
        }
        finally
        {
            try
            {
                if(bw !=null)
                    bw.close();
            }catch (Exception ex)
            {
                System.out.println("Error in closing the BufferedWriter"+ex);
            }
        }
    }
}
