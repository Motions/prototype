{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE LambdaCase #-}
import Bio.Motions.Prototype as Prototype
import Control.Monad.State.Strict
import Control.Monad.Except
import Options.Applicative
import System.IO
import System.Random
import Control.Lens
import Control.Monad.Random
import qualified Data.Serialize as Serialize
import qualified Data.ByteString as BS

data Settings
    = Initialize InitializeSettings
    | Simulate SimulateSettings

data InitializeSettings = InitializeSettings
    { settingsLaminsFile :: FilePath
    , settingsBindersFile :: FilePath
    , settingsRadius :: Int
    , settingsNumBinders :: Int
    , settingsChainLength :: Int
    , settingsInitOutputFile :: FilePath
    }

data SimulateSettings = SimulateSettings
    { settingsPDBFile :: FilePath
    , settingsInputFile :: FilePath
    , settingsOutputFile :: FilePath
    , settingsNumSteps :: Int
    , settingsWriteIntermediateStates :: Bool
    }

data Counter = Counter
    { _counterNumSteps :: !Int
    , _counterNumFrames :: !Int
    }
makeLenses ''Counter

initializeParser :: Parser InitializeSettings
initializeParser = InitializeSettings
    <$> strArgument
        (metavar "LAMIN-BSITES"
        <> help "File containing the lamin binding sites")
    <*> strArgument
        (metavar "BINDER-BSITES"
        <> help "File containing the regular binding sites")
    <*> option auto
        (long "radius"
        <> short 'r'
        <> metavar "RADIUS"
        <> help "Radius of the bounding sphere"
        <> value 10)
    <*> option auto
        (long "num-binders"
        <> short 'b'
        <> metavar "BINDERS"
        <> help "Number of binders"
        <> value 256)
    <*> option auto
        (long "chain-length"
        <> short 'l'
        <> help "Chain length"
        <> value 512)
    <*> strOption
        (long "output"
        <> short 'o'
        <> metavar "OUTPUT-FILE"
        <> help "Output file")

simulateParser :: Parser SimulateSettings
simulateParser = SimulateSettings
    <$> strOption
        (long "pdb"
        <> short 'p'
        <> metavar "PDB-FILE"
        <> help "PDB output file")
    <*> strOption
        (long "input"
        <> short 'i'
        <> metavar "INPUT-FILE"
        <> help "Input file")
    <*> strOption
        (long "output"
        <> short 'o'
        <> metavar "OUTPUT-FILE"
        <> help "Output file")
    <*> option auto
        (long "steps"
        <> short 's'
        <> metavar "STEPS"
        <> help "Number of simulaton steps"
        <> value 100000)
    <*> switch
        (long "intermediate-states"
        <> short 'm'
        <> help "Whether to write the intermediate states to the output file")

parser :: Parser Settings
parser = subparser
    $  command "init" (info (Initialize <$> initializeParser)
        (progDesc "Initialize the simulation"))
    <> command "run" (info (Simulate <$> simulateParser)
        (progDesc "Run the simulation"))

makeInput :: InitializeSettings -> IO Input
makeInput InitializeSettings{..} = do
    let inputChainLength = settingsChainLength
        inputRadius = settingsRadius
        inputNumBinders = settingsNumBinders
    inputLamins <- map read . lines <$> readFile settingsLaminsFile
    inputBinders <- map read . lines <$> readFile settingsBindersFile
    return Input{..}

runAndWrite :: (MonadSimulation m, MonadIO m) => Maybe Handle -> StateT Counter m ()
runAndWrite handle = do
    oldEnergy <- lift $ gets energy
    lift simulateStep
    newEnergy <- lift $ gets energy

    when (oldEnergy /= newEnergy) $ pushPDB handle

    counterNumSteps += 1

pushPDB :: (MonadSimulation m, MonadIO m) => Maybe Handle -> StateT Counter m ()
pushPDB Nothing = return ()
pushPDB (Just handle) = do
    st@SimulationState{..} <- lift get

    headerSequenceNumber <- use counterNumFrames
    headerStep <- use counterNumSteps
    let headerTitle = "chromosome;bonds=" ++ show energy

    liftIO $ do
        putStrLn $ "gyration radius: " ++ show gyrationRadius
        putStrLn $ "energy:          " ++ show energy
        writePDB handle Header{..} st
        hPutStrLn handle "END"

    counterNumFrames += 1

serialize :: FilePath -> SimulationState -> StdGen -> IO ()
serialize file st gen = BS.writeFile file $ Serialize.encode (st, show gen)

run :: Settings -> IO ()
run (Initialize settings) = do
    input <- makeInput settings
    (est, gen) <- newStdGen >>= runRandT (runExceptT $ initialize input)
    case est of
        Left err -> hPutStrLn stderr err
        Right st -> serialize (settingsInitOutputFile settings) st gen

run (Simulate SimulateSettings{..}) = withFile settingsPDBFile WriteMode $ \pdbFile ->
    Serialize.decode <$> BS.readFile settingsInputFile >>= \case
        Left err -> hPutStrLn stderr err
        Right (st, read -> gen) -> do
            (st', gen') <- flip runRandT gen $ flip execStateT st $ flip execStateT (Counter 0 0) $ do
                    pushPDB $ Just pdbFile
                    replicateM_ settingsNumSteps $ runAndWrite $
                        guard settingsWriteIntermediateStates >> Just pdbFile
            serialize settingsOutputFile st' gen'

main :: IO ()
main = run =<< execParser (info (helper <*> parser)
        (fullDesc
        <> progDesc "Perform a MCMC simulation of chromatine movements"))
