#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "data_management/output_connector_data.h"
#include "output_connector_data_transfer/output_connection.h"
#include "output_writer/generic.h"
#include "utility/string_utility.h"

//! iterate over all nested solvers
void SolverStructureVisualizer::
generateDiagram(std::stringstream &result, std::vector<std::tuple<int,int,int,SolverStructureVisualizer::solver_t::OutputConnection::output_connection_t>> &connectionLines,
      std::vector<int> &slotLineNos, int depth, bool isFirstChildSolver, bool isLastChildSolver)
{
  if (!currentSolver_)
    return;

  if (!enabled_)
    return;

  slotLineNos.clear();

  std::stringstream lineStart;
  for (int i = 0; i < depth; i++)
    lineStart << "│ ";

  result << lineStart.str();
  if (isLastChildSolver && currentSolver_->outputSlots.empty() && currentSolver_->outputConnections.empty())
    result << "└";
  else
    result << "├";
    //result << "┌";
  result << "── " << currentSolver_->name << "\n";

  if (currentSolver_->description != "")
  {
    result << lineStart.str() << "│   (\"" << currentSolver_->description << "\")" << "\n";
  }

  // print the output slots of the solver, store the line nos of each slot to slotLineNos
  if (!currentSolver_->outputSlots.empty())
  {
    result << lineStart.str() << "│  output slots: \n";

    for (int i = 0; i < currentSolver_->outputSlots.size(); i++)
    {
      std::stringstream s;
      s << lineStart.str() << "│  " << currentSolver_->outputSlots[i].fieldVariableName << "."
        << currentSolver_->outputSlots[i].componentName;

      std::string outputSlotString = s.str();

      int currentLineLength = StringUtility::stringLength(outputSlotString)+2;

      const int requiredLineLength = 45;

      s << std::string(requiredLineLength - currentLineLength, ' ')
        << "── ¤" << i << "\n";
      result << s.str();

      // store line no
      std::string currentResult = result.str();
      int currentLineNo = std::count(currentResult.begin(), currentResult.end(), '\n')-1;

      //LOG(DEBUG) << "current Result [" << currentResult << "] contains " << currentLineNo << " line breaks";
      slotLineNos.push_back(currentLineNo);
    }
    result << lineStart.str() << "│\n";
  }

  // print the output slot connections of the solver
  std::vector<solver_t::OutputConnection> &outputConnections = currentSolver_->outputConnections;
  if (!outputConnections.empty())
  {
    result << lineStart.str() << "│  slot connections: \n";
    for (int i = 0; i < outputConnections.size(); i++)
    {
      if (outputConnections[i].type == solver_t::OutputConnection::output_connection_t::ba)
      {
        result << lineStart.str() << "│  " <<  outputConnections[i].toSlot << "¤ <- ¤" << outputConnections[i].fromSlot << "\n";
      }
      else
      {
        result << lineStart.str() << "│  " <<  outputConnections[i].fromSlot << "¤";

        switch (outputConnections[i].type)
        {
        case solver_t::OutputConnection::output_connection_t::ba:
        case solver_t::OutputConnection::output_connection_t::ab:
          result << " -> ";
          break;
        case solver_t::OutputConnection::output_connection_t::bidirectionalReuse:
          result << " <=> ";
          break;
        case solver_t::OutputConnection::output_connection_t::bidirectionalCopy:
          result << " <-> ";
          break;

        };

        result << "¤" << outputConnections[i].toSlot << "\n";
      }
    }
    result << lineStart.str() << "│\n";
  }

  // get the line nos of the slots for nested solvers 0 and 1
  std::vector<std::shared_ptr<solver_t>> &children = currentSolver_->children;
  std::vector<std::vector<int>> slotLineNosTerm(children.size());
  int i = 0;
  // loop over child solvers
  for (std::vector<std::shared_ptr<solver_t>>::iterator iter = children.begin();
       iter != children.end(); iter++, i++)
  {
    currentSolver_ = (*iter);

    bool isFirstChildSolver = i == 0;
    bool isLastChildSolver = i == children.size()-1;

    // recursively call generateDiagram on subsolvers
    generateDiagram(result, connectionLines, slotLineNosTerm[i], depth+1, isFirstChildSolver, isLastChildSolver);

    // add end of subsolver
    if (isLastChildSolver)
      result << lineStart.str() << "└\n";
    else
      result << lineStart.str() << "│\n";
  }

  // add connections between slots of nested solvers to connectionLines
  for (int i = 0; i < outputConnections.size(); i++)
  {
    if (slotLineNosTerm.size() >= 2)
    {
      if ((outputConnections[i].type == solver_t::OutputConnection::output_connection_t::ab
        || outputConnections[i].type == solver_t::OutputConnection::output_connection_t::bidirectionalCopy
        || outputConnections[i].type == solver_t::OutputConnection::output_connection_t::bidirectionalReuse)
        && slotLineNosTerm[0].size() > outputConnections[i].fromSlot
        && slotLineNosTerm[1].size() > outputConnections[i].toSlot
      )
      {
        int lineNoFrom = slotLineNosTerm[0][outputConnections[i].fromSlot];
        int lineNoTo = slotLineNosTerm[1][outputConnections[i].toSlot];
        solver_t::OutputConnection::output_connection_t connectionType = outputConnections[i].type;

        connectionLines.push_back(std::make_tuple(
          lineNoFrom, lineNoTo, 0, connectionType
        ));
      }
      else if (outputConnections[i].type == solver_t::OutputConnection::output_connection_t::ba
        && slotLineNosTerm[1].size() > outputConnections[i].fromSlot
        && slotLineNosTerm[0].size() > outputConnections[i].toSlot
      )
      {
        int lineNoFrom = slotLineNosTerm[0][outputConnections[i].toSlot];
        int lineNoTo = slotLineNosTerm[1][outputConnections[i].fromSlot];
        solver_t::OutputConnection::output_connection_t connectionType = outputConnections[i].type;

        connectionLines.push_back(std::make_tuple(
          lineNoFrom, lineNoTo, 0, connectionType
        ));
      }
      else
      {
        LOG(DEBUG) << "slotLineNosTerm is not good. fromSlot: " << outputConnections[i].fromSlot
          << ", toSlot: " << outputConnections[i].toSlot << ", slot line nos term0: " << slotLineNosTerm[0]
          << ", term1: " << slotLineNosTerm[1];
      }
    }
    else
    {
      LOG(DEBUG) << "slotLineNosTerm is too short: " << slotLineNosTerm.size() << ", but there are " << outputConnections.size() << " outputConnections.";
    }
  }

  // if the own solver does not have output slots but it has only a single child, use their output slots,
  // this is true for MultipleInstances
  if (slotLineNos.empty() && children.size() == 1)
  {
    slotLineNos = slotLineNosTerm[0];
  }
}

//! produce the resulting file
void SolverStructureVisualizer::
writeDiagramFile(std::string filename)
{
  LOG(DEBUG) << "SolverStructureVisualizer::writeDiagramFile() nDisableCalls_: " << nDisableCalls_ << ", enabled: " << enabled_ << ", currently at \"" << currentSolver_->name << "\".";

  if (!enabled_)
    return;

  // only produce file on rank 0
  if (DihuContext::ownRankNoCommWorld() != 0)
    return;


  // collect all information

  // backup currentSolver pointer
  std::shared_ptr<solver_t> currentSolverBeforePrint = currentSolver_;

  // set to root of nested solvers to start printing from there
  currentSolver_ = solverRoot_;

  //! connection lines contains the following information <lineNoFrom, lineNoTo, lineColumn, lineType>
  // lineColumn is the horizontal position of the vertical data connection line
  std::vector<std::tuple<int,int,int,SolverStructureVisualizer::solver_t::OutputConnection::output_connection_t>> connectionLines;
  std::vector<int> slotLineNos;
  std::stringstream diagramWithoutDataLines;

  // call generateDiagram on the nested solvers
  generateDiagram(diagramWithoutDataLines, connectionLines, slotLineNos);

  // restore currentSolver
  currentSolver_ = currentSolverBeforePrint;

  // set the "lineColumn" entry in connectionLines
  // loop over connectionLines
  for (int i = 0; i < connectionLines.size(); i++)
  {
    int startLineNo = std::get<0>(connectionLines[i]);
    int endLineNo = std::get<1>(connectionLines[i]);
    int lineColumn = 0;

    // count number of starting lines inside this line
    for (int j = 0; j < connectionLines.size(); j++)
    {
      if (i == j)
        continue;
      int startLineNo2 = std::get<0>(connectionLines[j]);
      int endLineNo2 = std::get<1>(connectionLines[j]);

      if (startLineNo < startLineNo2 && startLineNo2 < endLineNo)
        lineColumn++;
      else if (startLineNo == startLineNo2 && endLineNo2 < endLineNo)
        lineColumn++;
      else if (startLineNo < startLineNo2 && endLineNo2 == endLineNo)
        lineColumn++;
    }
    std::get<2>(connectionLines[i]) = lineColumn;
  }

  LOG(DEBUG) << "connectionLines: " << connectionLines;

  std::size_t pos = 0;
  std::string diagramWithoutDataLinesString = diagramWithoutDataLines.str();
  std::string diagram;
  diagram += "Solver structure: \n\n";

  // print actual diagram, into diagram
  // loop over lines
  for (int currentLineNo = 0; pos < diagramWithoutDataLinesString.length(); currentLineNo++)
  {
    // determine current line by iterating over substring between "\n"s
    int nextPos = diagramWithoutDataLinesString.find("\n",pos+1);
    std::string line;

    if (pos == 0)
      pos = -1;

    if (nextPos != std::string::npos)
    {
      line = diagramWithoutDataLinesString.substr(pos+1, nextPos-pos-1);
    }
    else
    {
      line = diagramWithoutDataLinesString.substr(pos+1);
    }

    pos = nextPos;

    // here, line contains the current line that should be handled without "\n"

    // find out vertial data slot connection lines that go through this output line
    std::vector<std::tuple<int,bool,bool,solver_t::OutputConnection::output_connection_t>> verticalLinesInCurrentRow;    //! <lineColumn,thisLineStartsHere,thisLineEndsHere,type>, lineColumn is the column in rhe current row where the vertical lines passes through
    bool lineStartsHere = false;
    bool lineEndsHere = false;
    solver_t::OutputConnection::output_connection_t lineType = solver_t::OutputConnection::output_connection_t::ab;
    int lineColumn = 0;
#if 0
    // get maximum lineColumn
    int maximumLinePosition = 0;
    for (int i = 0; i < connectionLines.size(); i++)
    {
      maximumLinePosition = std::max((int)(std::get<2>(connectionLines[i])), maximumLinePosition);
    }
#endif
    // loop over all connectionLines and see if they pass through the current line
    // collect all those rows into data structure verticalLinesInCurrentRow
    for (int i = 0; i < connectionLines.size(); i++)
    {
      // start and end row of current connection line
      int startLineNo = std::get<0>(connectionLines[i]);
      int endLineNo = std::get<1>(connectionLines[i]);

      // if the connection line starts or ends in the current row, there will be a horizontal line
      if (currentLineNo == startLineNo)
      {
        lineType = std::get<3>(connectionLines[i]);
        lineStartsHere = true;
        lineColumn = std::get<2>(connectionLines[i]);
        verticalLinesInCurrentRow.push_back(std::make_tuple(lineColumn, true, false, lineType));
      }
      else if (currentLineNo == endLineNo)
      {
        lineType = std::get<3>(connectionLines[i]);
        lineEndsHere = true;
        lineColumn = std::get<2>(connectionLines[i]);
        verticalLinesInCurrentRow.push_back(std::make_tuple(lineColumn, false, true, lineType));
      }
      else if (startLineNo < currentLineNo && currentLineNo < endLineNo)
      {
        // vertical line passes through current row
        solver_t::OutputConnection::output_connection_t lineType = std::get<3>(connectionLines[i]);
        int lineColumn = std::get<2>(connectionLines[i]);
        verticalLinesInCurrentRow.push_back(std::make_tuple(lineColumn, false, false, lineType));
      }
    }

    // sort all found lines according to their column
    std::sort(verticalLinesInCurrentRow.begin(), verticalLinesInCurrentRow.end(), [](
      const std::tuple<int,bool,bool,solver_t::OutputConnection::output_connection_t> &a,
      const std::tuple<int,bool,bool,solver_t::OutputConnection::output_connection_t> &b)
    {
      return std::get<0>(a) < std::get<0>(b);
    });

    VLOG(1) << "currentLineNo " << currentLineNo << ", line=[" << line << "] has vertical lines: " << verticalLinesInCurrentRow;

    // fill current line with spaces such that it has the correct length
    const int requiredLineLength = 49;
    int lineLength = StringUtility::stringLength(line);

    // determine fill character for horizontal line
    std::string fillCharacter = " ";
    if (lineStartsHere || lineEndsHere)
    {
      if (lineType == solver_t::OutputConnection::output_connection_t::bidirectionalReuse || lineType == solver_t::OutputConnection::output_connection_t::bidirectionalCopy)
      {
        fillCharacter = "═";
      }
      else
      {
        fillCharacter = "─";
      }
    }

    // if the line is yet too short, add fillCharacter until is has the correct length
    if (requiredLineLength > lineLength)
    {
      for (int i = 0; i < requiredLineLength - lineLength; i++)
      {
        // starting tip of arrow if the line starts or ends in the current row and has the according type
        if (i == 0 && lineStartsHere && lineType == solver_t::OutputConnection::output_connection_t::ba)
        {
          line += "<";
        }
        else if (i == 0 && lineEndsHere && lineType == solver_t::OutputConnection::output_connection_t::ab)
        {
          line += "<";
        }
        else
        {
          line += fillCharacter;
        }
      }
    }

    // add the code line to the diagram, so far without the vertical lines
    diagram += line;

    // loop over the vertical lines in the current row
    int linesIndex = 0;
    if (!verticalLinesInCurrentRow.empty())
    {
      // loop over columns in the current row, current column is requiredLineLength + 2*i, i.e. it advances always by 2 characters
      for (int i = 0; i <= std::get<0>(verticalLinesInCurrentRow.back()); i++)
      {

        // determine events where line starts or ends
        bool lineCrosses = (lineStartsHere || lineEndsHere) && i < lineColumn;
        bool lineCurveStart = lineStartsHere && i == lineColumn;
        bool lineCurveEnd = lineEndsHere && i == lineColumn;

        // determine space character which is a horizontal line if the connection starts or ends here
        std::string space = " ";
        if (lineCrosses || ( (lineStartsHere || lineEndsHere) && i <= lineColumn))
        {
          if (lineType == solver_t::OutputConnection::output_connection_t::bidirectionalReuse
              || lineType == solver_t::OutputConnection::output_connection_t::bidirectionalCopy)
          {
            space = "═";
          }
          else
          {
            space = "─";
          }
        }

        // if there is a vertical line passing at the current position
        if (i == std::get<0>(verticalLinesInCurrentRow[linesIndex]))
        {
          // depending on the type of the vertical line choose the right sign
          switch (std::get<3>(verticalLinesInCurrentRow[linesIndex]))
          {
          case solver_t::OutputConnection::output_connection_t::bidirectionalReuse:
            diagram += space;
            if (lineCrosses && !std::get<1>(verticalLinesInCurrentRow[linesIndex]) && !std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              diagram += "╬";
            }
            else if (lineCurveStart)
            {
              diagram += "╗";
            }
            else if (lineCurveEnd)
            {
              diagram += "╝";
            }
            else if (std::get<1>(verticalLinesInCurrentRow[linesIndex]))
            {
              diagram += "╦";
            }
            else if (std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              diagram += "╩";
            }
            else
            {
              diagram += "║";
            }
            break;
          case solver_t::OutputConnection::output_connection_t::ab:
          case solver_t::OutputConnection::output_connection_t::ba:
          case solver_t::OutputConnection::output_connection_t::bidirectionalCopy:
            diagram += space;
            if (lineCrosses && !std::get<1>(verticalLinesInCurrentRow[linesIndex]) && !std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              diagram += "┼";
            }
            else if (lineCurveStart)
            {
              diagram += "┐";
            }
            else if (lineCurveEnd)
            {
              diagram += "┘";
            }
            else if (std::get<1>(verticalLinesInCurrentRow[linesIndex]))
            {
              diagram += "┬";
            }
            else if (std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              diagram += "┴";
            }
            else
            {
              diagram += "│";
            }
            break;
          };
          linesIndex++;
        }
        else
        {
          // there is no vertical line passing through the current position, add two spaces
          diagram += space;
          diagram += space;
        }
      }
    }

    // add endline at the end of the current line
    diagram += "\n";
  }

  diagram += "connection types:\n";
  diagram += "  ═══ ... reuse field variable, no copy\n";
  diagram += "  ──> ... copy data in direction of arrow\n";

  //LOG(DEBUG) << "solver structure without connections:\n" << diagramWithoutDataLines.str() << "\n";
  LOG(DEBUG) << "solver structure with connections:\n" << diagram << "\n";

  // output diagram to a text file
  std::ofstream file;
  OutputWriter::Generic::openFile(file, filename);

  if (!file.is_open())
  {
    LOG(FATAL) << "Could not write to file \"" << filename << "\".";
  }

  // output diagram text
  file << diagram;
  file.close();

  LOG(INFO) << "File \"" << filename << "\" written.";
}
