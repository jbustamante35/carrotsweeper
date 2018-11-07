import org.json.simple.JSONArray
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
midline = JSONArray;
%%
midline.put(num2str(rand(100,2)));
what = JSONObject();
what.put('hello',midline);
what.toString()
%%
file = File('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Arabidopsis/arabidopsis_auto/midline.json');
file.createNewFile();
FileWriter(file);

fileWriter.write(what.toJSONString());
fileWriter.flush();
fileWriter.close();